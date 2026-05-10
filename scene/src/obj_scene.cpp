#include "scene/obj_scene.hpp"

#include <cmath>
#include <fstream>
#include <nlohmann/json.hpp>
#include <sstream>
#include <stdexcept>

#include "spdlog/spdlog.h"
#include "stb_image.h"

WOS_NAMESPACE_OPEN_SCOPE

namespace {

template<typename ScalarType>
ScalarType
extractChannel(const Vector<3>& val, ChannelMode mode)
{
    if constexpr (std::is_same_v<ScalarType, Scalar<3>>) {
        return Scalar<3>(val[0], val[1], val[2]);
    } else {
        switch (mode) {
            case ChannelMode::Red:
                return ScalarType(val[0]);
            case ChannelMode::Green:
                return ScalarType(val[1]);
            case ChannelMode::Blue:
                return ScalarType(val[2]);
            default:
                return ScalarType(val[0]);
        }
    }
}

template<typename ScalarType>
double
scalarL1Norm(const ScalarType& s)
{
    if constexpr (std::is_same_v<ScalarType, Scalar<1>>)
        return std::abs(s.value());
    else
        return std::abs(s[0]) + std::abs(s[1]) + std::abs(s[2]);
}

template<typename ScalarType>
ScalarType
scalarNegate(const ScalarType& v)
{
    if constexpr (std::is_same_v<ScalarType, Scalar<1>>)
        return ScalarType(-v.value());
    else
        return Scalar<3>(-v[0], -v[1], -v[2]);
}

template<typename ScalarType>
ScalarType
scalarLerp(const ScalarType& a, const ScalarType& b, double t)
{
    if constexpr (std::is_same_v<ScalarType, Scalar<1>>)
        return ScalarType((1.0 - t) * a.value() + t * b.value());
    else
        return Scalar<3>((1.0 - t) * a[0] + t * b[0], (1.0 - t) * a[1] + t * b[1], (1.0 - t) * a[2] + t * b[2]);
}

} // namespace

// ============================================================================
// Constructor / Destructor
// ============================================================================

template<typename ScalarType, int DIM>
ObjScene<ScalarType, DIM>::ObjScene(const std::string& sceneDir, ChannelMode channel)
  : channel_(channel)
{
    this->interest_scene_ = std::make_unique<fcpw::Scene<DIM>>();

    loadConfig(sceneDir + "/config.json");
    if constexpr (DIM == 3)
        loadBoundaryTextures(sceneDir);
    loadGeometryAndBoundary(sceneDir + "/mesh.obj", sceneDir + "/boundary.txt");
    loadInterest(sceneDir + "/interest.obj");
    if constexpr (DIM == 2)
        loadSource(sceneDir + "/source.bmp");
}

template<typename ScalarType, int DIM>
ObjScene<ScalarType, DIM>::~ObjScene()
{
    if constexpr (DIM == 2) {
        if (texData_)
            stbi_image_free(texData_);
    } else {
        if (dTexData_)
            stbi_image_free(dTexData_);
        if (nTexData_)
            stbi_image_free(nTexData_);
    }
}

// ============================================================================
// loadConfig
// ============================================================================

template<typename ScalarType, int DIM>
void
ObjScene<ScalarType, DIM>::loadConfig(const std::string& path)
{
    std::ifstream f(path);
    if (!f.is_open())
        return;

    nlohmann::json j;
    f >> j;

    if (j.contains("dimension") && j["dimension"] != DIM) {
        spdlog::error("Scene dimension mismatch: expected {}, got {}", DIM, j["dimension"].get<int>());
        throw std::runtime_error("Scene dimension mismatch");
    }

    if (j.contains("interior_only")) {
        this->interior_only_ = j["interior_only"];
    }

    if (j.contains("render_range")) {
        for (int i = 0; i < DIM; ++i) {
            this->renderMin_[i] = j["render_range"]["min"][i];
            this->renderMax_[i] = j["render_range"]["max"][i];
        }
    }
    if (j.contains("source_range")) {
        for (int i = 0; i < DIM; ++i) {
            sourceMin_[i] = j["source_range"]["min"][i];
            sourceMax_[i] = j["source_range"]["max"][i];
        }
    }
    if (j.contains("source_intensity"))
        sourceIntensity_ = j["source_intensity"];
    if (j.value("is_screened", false)) {
        this->is_screened_ = true;
        if (j.contains("screened_kappa"))
            this->screened_kappa_ = j["screened_kappa"];
    }

    if (j.contains("solver_plane")) {
        auto& sp = j["solver_plane"];
        for (int i = 0; i < 3; ++i) {
            solverPlaneOrigin_[i] = sp["origin"][i];
            solverPlaneUAxis_[i] = sp["u_axis"][i];
            solverPlaneVAxis_[i] = sp["v_axis"][i];
        }
        solverPlaneWidth_ = sp.value("width", 2.0);
        solverPlaneHeight_ = sp.value("height", 2.0);

        if (!j.contains("render_range")) {
            double hw = 0.5 * solverPlaneWidth_;
            double hh = 0.5 * solverPlaneHeight_;
            for (int d = 0; d < DIM; ++d) {
                double du = hw * std::abs(solverPlaneUAxis_[d]);
                double dv = hh * std::abs(solverPlaneVAxis_[d]);
                this->renderMin_[d] = solverPlaneOrigin_[d] - du - dv;
                this->renderMax_[d] = solverPlaneOrigin_[d] + du + dv;
            }
        }
    }
}

// ============================================================================
// loadGeometryAndBoundary — inlined 2D / 3D
// ============================================================================

template<typename ScalarType, int DIM>
void
ObjScene<ScalarType, DIM>::loadGeometryAndBoundary(const std::string& objPath, const std::string& boundaryPath)
{
    if constexpr (DIM == 2) {
        // ---- 2D: line segments from OBJ + boundary.txt ----
        using Vec2 = fcpw::Vector<2>;

        std::vector<Vec2> rawVertices;
        std::vector<Vec2> dVertices, nVertices;
        std::vector<Eigen::Vector2i> dEdges, nEdges, allEdges;
        std::vector<int> dIndex, nIndex;

        {
            std::ifstream f(objPath);
            if (!f.is_open()) {
                spdlog::error("Failed to open OBJ: {}", objPath);
                throw std::runtime_error("Failed to open OBJ: " + objPath);
            }
            std::string line;
            while (std::getline(f, line)) {
                if (line.rfind("v ", 0) != 0)
                    continue;
                std::istringstream iss(line.substr(2));
                float x, y, z;
                iss >> x >> y >> z;
                Vec2 p(x, -z);
                rawVertices.push_back(p);
                dIndex.push_back(-1);
                nIndex.push_back(-1);
            }
        }

        {
            std::ifstream bf(boundaryPath);
            if (!bf.is_open()) {
                spdlog::error("Failed to open boundary file: {}", boundaryPath);
                throw std::runtime_error("Failed to open boundary file: " + boundaryPath);
            }
            std::string line;
            while (std::getline(bf, line)) {
                std::istringstream iss(line);
                int type, v[2];
                if (!(iss >> v[0] >> v[1] >> type))
                    continue;
                allEdges.emplace_back(v[0], v[1]);

                double d_s[3], d_e[3], n_s[3], n_e[3];
                iss >> d_s[0] >> d_s[1] >> d_s[2] >> d_e[0] >> d_e[1] >> d_e[2] >> n_s[0] >> n_s[1] >> n_s[2] >>
                  n_e[0] >> n_e[1] >> n_e[2];

                EdgeBoundaryData data;
                int isDouble = 0;
                if (iss >> isDouble) {
                    data.isDoubleSided = (isDouble != 0);
                    if (data.isDoubleSided) {
                        double d_s2[3], d_e2[3], n_s2[3], n_e2[3];
                        iss >> d_s2[0] >> d_s2[1] >> d_s2[2] >> d_e2[0] >> d_e2[1] >> d_e2[2] >> n_s2[0] >> n_s2[1] >>
                          n_s2[2] >> n_e2[0] >> n_e2[1] >> n_e2[2];
                        data.startVal2 = extractChannel<ScalarType>(Vector<3>(d_s2[0], d_s2[1], d_s2[2]), channel_);
                        data.endVal2 = extractChannel<ScalarType>(Vector<3>(d_e2[0], d_e2[1], d_e2[2]), channel_);
                    }
                }

                auto& targetData = (type == 0) ? dirichletData_ : neumannData_;
                auto& targetIndex = (type == 0) ? dIndex : nIndex;
                auto& targetVerts = (type == 0) ? dVertices : nVertices;
                auto& targetEdges = (type == 0) ? dEdges : nEdges;
                const double* s = (type == 0) ? d_s : n_s;
                const double* e = (type == 0) ? d_e : n_e;

                if (type == 2) {
                    spdlog::error("ObjScene<2>: Mixed Dirichlet/Neumann boundary not supported.");
                    throw std::runtime_error("ObjScene<2>: Mixed Dirichlet/Neumann boundary not supported.");
                }

                for (int i = 0; i < 2; ++i) {
                    if (targetIndex[v[i]] == -1) {
                        targetIndex[v[i]] = static_cast<int>(targetVerts.size());
                        targetVerts.push_back(rawVertices[v[i]]);
                    }
                }
                targetEdges.emplace_back(targetIndex[v[0]], targetIndex[v[1]]);
                data.startVal = extractChannel<ScalarType>(Vector<3>(s[0], s[1], s[2]), channel_);
                data.endVal = extractChannel<ScalarType>(Vector<3>(e[0], e[1], e[2]), channel_);
                targetData.push_back(data);
            }
        }

        if (!dEdges.empty()) {
            this->dirichlet_scene_ = std::make_unique<fcpw::Scene<2>>();
            this->dirichlet_scene_->setObjectCount(1);
            this->dirichlet_scene_->setObjectVertices(dVertices, 0);
            this->dirichlet_scene_->setObjectLineSegments(dEdges, 0);
            this->dirichlet_scene_->build(fcpw::AggregateType::Bvh_SurfaceArea, true);
            active_scenes_.push_back(this->dirichlet_scene_.get());
        }
        if (!nEdges.empty()) {
            this->neumann_scene_ = std::make_unique<fcpw::Scene<2>>();
            this->neumann_scene_->setObjectCount(1);
            this->neumann_scene_->setObjectVertices(nVertices, 0);
            this->neumann_scene_->setObjectLineSegments(nEdges, 0);
            this->neumann_scene_->computeSilhouettes([](float, int) { return false; });
            this->neumann_scene_->build(fcpw::AggregateType::Bvh_SurfaceArea, true);
            active_scenes_.push_back(this->neumann_scene_.get());
        }
        if (!allEdges.empty()) {
            this->full_scene_ = std::make_unique<fcpw::Scene<2>>();
            this->full_scene_->setObjectCount(1);
            this->full_scene_->setObjectVertices(rawVertices, 0);
            this->full_scene_->setObjectLineSegments(allEdges, 0);
            this->full_scene_->build(fcpw::AggregateType::Bvh_SurfaceArea, true);
        }
    } else {
        // ---- 3D: triangles from OBJ + boundary.txt ----
        using Vec3 = fcpw::Vector<3>;

        std::vector<int> faceTypes;
        {
            std::ifstream f(boundaryPath);
            if (f.is_open()) {
                int fidx, type;
                while (f >> fidx >> type) {
                    if (fidx >= static_cast<int>(faceTypes.size()))
                        faceTypes.resize(fidx + 1, 0);
                    faceTypes[fidx] = type;
                }
            }
        }

        std::vector<Vec3> rawVertices;
        std::vector<Vector<2>> rawUVs;
        std::vector<Vec3> dVertices, nVertices;
        std::vector<Eigen::Vector3i> dTris, nTris, allTris;
        std::vector<int> dIndex, nIndex;

        std::ifstream f(objPath);
        if (!f.is_open()) {
            spdlog::error("Failed to open OBJ: {}", objPath);
            throw std::runtime_error("Failed to open OBJ: " + objPath);
        }

        std::string line;
        int faceIdx = 0;
        while (std::getline(f, line)) {
            if (line.rfind("v ", 0) == 0) {
                std::istringstream iss(line.substr(2));
                float x, y, z;
                iss >> x >> y >> z;
                Vec3 p(x, -z, y);
                rawVertices.push_back(p);
                dIndex.push_back(-1);
                nIndex.push_back(-1);
            } else if (line.rfind("vt ", 0) == 0) {
                std::istringstream iss(line.substr(3));
                float u, v;
                iss >> u >> v;
                rawUVs.emplace_back(u, v);
            } else if (line.rfind("f ", 0) == 0) {
                std::istringstream iss(line.substr(2));
                std::string part;
                struct VertexIdx
                {
                    int v = -1, vt = -1;
                };
                std::vector<VertexIdx> faceVerts;
                while (iss >> part) {
                    VertexIdx vi;
                    size_t s1 = part.find('/');
                    if (s1 == std::string::npos) {
                        vi.v = std::stoi(part) - 1;
                    } else {
                        vi.v = std::stoi(part.substr(0, s1)) - 1;
                        size_t s2 = part.find('/', s1 + 1);
                        std::string vtStr =
                          (s2 == std::string::npos) ? part.substr(s1 + 1) : part.substr(s1 + 1, s2 - (s1 + 1));
                        if (!vtStr.empty())
                            vi.vt = std::stoi(vtStr) - 1;
                    }
                    faceVerts.push_back(vi);
                }

                int type = (faceIdx < static_cast<int>(faceTypes.size())) ? faceTypes[faceIdx] : 0;
                faceIdx++;

                if (faceVerts.size() < 3)
                    continue;

                for (size_t i = 1; i + 1 < faceVerts.size(); ++i) {
                    int v[3] = { faceVerts[0].v, faceVerts[i].v, faceVerts[i + 1].v };
                    allTris.emplace_back(v[0], v[1], v[2]);

                    int vt[3] = { faceVerts[0].vt, faceVerts[i].vt, faceVerts[i + 1].vt };
                    Vector<2> uv[3];
                    for (int k = 0; k < 3; ++k)
                        uv[k] =
                          (vt[k] >= 0 && vt[k] < static_cast<int>(rawUVs.size())) ? rawUVs[vt[k]] : Vector<2>(0, 0);

                    if (type == 0) {
                        for (int k = 0; k < 3; ++k) {
                            if (dIndex[v[k]] == -1) {
                                dIndex[v[k]] = static_cast<int>(dVertices.size());
                                dVertices.push_back(rawVertices[v[k]]);
                            }
                        }
                        dTris.emplace_back(dIndex[v[0]], dIndex[v[1]], dIndex[v[2]]);
                        dirichletUVs_.push_back({ uv[0], uv[1], uv[2] });
                    } else if (type == 1) {
                        for (int k = 0; k < 3; ++k) {
                            if (nIndex[v[k]] == -1) {
                                nIndex[v[k]] = static_cast<int>(nVertices.size());
                                nVertices.push_back(rawVertices[v[k]]);
                            }
                        }
                        nTris.emplace_back(nIndex[v[0]], nIndex[v[1]], nIndex[v[2]]);
                        neumannUVs_.push_back({ uv[0], uv[1], uv[2] });
                    } else if (type == 2) {
                        spdlog::error("ObjScene<3>: Mixed Dirichlet/Neumann boundary not supported.");
                        throw std::runtime_error("ObjScene<3>: Mixed Dirichlet/Neumann boundary not supported.");
                    }
                }
            }
        }

        if (!dTris.empty()) {
            this->dirichlet_scene_ = std::make_unique<fcpw::Scene<3>>();
            this->dirichlet_scene_->setObjectCount(1);
            this->dirichlet_scene_->setObjectVertices(dVertices, 0);
            this->dirichlet_scene_->setObjectTriangles(dTris, 0);
            this->dirichlet_scene_->build(fcpw::AggregateType::Bvh_SurfaceArea, true);
            active_scenes_.push_back(this->dirichlet_scene_.get());
        }
        if (!nTris.empty()) {
            this->neumann_scene_ = std::make_unique<fcpw::Scene<3>>();
            this->neumann_scene_->setObjectCount(1);
            this->neumann_scene_->setObjectVertices(nVertices, 0);
            this->neumann_scene_->setObjectTriangles(nTris, 0);
            this->neumann_scene_->computeSilhouettes([](float, int) { return false; });
            this->neumann_scene_->build(fcpw::AggregateType::Bvh_SurfaceArea, true, true);
            active_scenes_.push_back(this->neumann_scene_.get());
        }
        if (!allTris.empty()) {
            this->full_scene_ = std::make_unique<fcpw::Scene<3>>();
            this->full_scene_->setObjectCount(1);
            this->full_scene_->setObjectVertices(rawVertices, 0);
            this->full_scene_->setObjectTriangles(allTris, 0);
            this->full_scene_->build(fcpw::AggregateType::Bvh_SurfaceArea, true);
        }
    }
}

// ============================================================================
// loadInterest
// ============================================================================

template<typename ScalarType, int DIM>
void
ObjScene<ScalarType, DIM>::loadInterest(const std::string& path)
{
    if constexpr (DIM == 3) {
        this->has_interest_ = false;
    } else {
        std::ifstream f(path);
        if (!f.is_open()) {
            this->has_interest_ = false;
            return;
        }

        std::vector<fcpw::Vector<2>> vertices;
        std::vector<Eigen::Vector2i> segments;
        std::string line;
        while (std::getline(f, line)) {
            if (line.empty())
                continue;
            std::istringstream iss(line);
            std::string type;
            iss >> type;
            if (type == "v") {
                float x, y, z;
                iss >> x >> y >> z;
                vertices.emplace_back(x, -z);
            } else if (type == "l") {
                int v1, v2;
                iss >> v1 >> v2;
                segments.emplace_back(v1 - 1, v2 - 1);
            }
        }
        if (!segments.empty()) {
            this->interest_scene_->setObjectCount(1);
            this->interest_scene_->setObjectVertices(vertices, 0);
            this->interest_scene_->setObjectLineSegments(segments, 0);
            this->interest_scene_->build(fcpw::AggregateType::Bvh_SurfaceArea, true);
            this->has_interest_ = true;
        }
    }
}

// ============================================================================
// loadSource (2D only)
// ============================================================================

template<typename ScalarType, int DIM>
void
ObjScene<ScalarType, DIM>::loadSource(const std::string& path)
{
    if constexpr (DIM == 3)
        return;
    texData_ = stbi_load(path.c_str(), &texWidth_, &texHeight_, &texChannels_, 0);
    if (!texData_)
        spdlog::warn("Failed to load source texture: {}", path);
}

// ============================================================================
// loadBoundaryTextures (3D only)
// ============================================================================

template<typename ScalarType, int DIM>
void
ObjScene<ScalarType, DIM>::loadBoundaryTextures(const std::string& dir)
{
    if constexpr (DIM == 2)
        return;
    dTexData_ = stbi_load((dir + "/dirichlet.bmp").c_str(), &dTexWidth_, &dTexHeight_, &dTexChannels_, 0);
    nTexData_ = stbi_load((dir + "/neumann.bmp").c_str(), &nTexWidth_, &nTexHeight_, &nTexChannels_, 0);
}

// ============================================================================
// source (1-param override + 2-param implementation)
// ============================================================================

template<typename ScalarType, int DIM>
ScalarType
ObjScene<ScalarType, DIM>::source(const VectorType& p) const
{
    return source(p, nullptr);
}

template<typename ScalarType, int DIM>
ScalarType
ObjScene<ScalarType, DIM>::source(const VectorType& p, bool* is_zero) const
{
    if constexpr (DIM == 3) {
        if (is_zero)
            *is_zero = true;
        return ScalarType(0.0);
    }
    return sampleTexture(p, is_zero);
}

// ============================================================================
// sampleTexture (2D)
// ============================================================================

template<typename ScalarType, int DIM>
ScalarType
ObjScene<ScalarType, DIM>::sampleTexture(const VectorType& p, bool* is_zero) const
{
    if (is_zero)
        *is_zero = true;
    if (!texData_)
        return ScalarType(0.0);

    double u = (p[0] - sourceMin_[0]) / (sourceMax_[0] - sourceMin_[0]);
    double v = (p[1] - sourceMin_[1]) / (sourceMax_[1] - sourceMin_[1]);
    if (u < 0 || u > 1 || v < 0 || v > 1)
        return ScalarType(0.0);

    v = 1.0 - v;
    double x = u * (texWidth_ - 1), y = v * (texHeight_ - 1);
    int x0 = static_cast<int>(std::floor(x)), y0 = static_cast<int>(std::floor(y));
    int x1 = std::min(x0 + 1, texWidth_ - 1), y1 = std::min(y0 + 1, texHeight_ - 1);
    double dx = x - x0, dy = y - y0;

    auto getPixel = [&](int px, int py) -> Vector<3> {
        int idx = (py * texWidth_ + px) * texChannels_;
        double r = texData_[idx] / 255.0;
        double g = (texChannels_ > 1) ? texData_[idx + 1] / 255.0 : r;
        double b = (texChannels_ > 2) ? texData_[idx + 2] / 255.0 : 0.0;
        return Vector<3>(r, g, b);
    };

    Vector<3> v00 = getPixel(x0, y0), v10 = getPixel(x1, y0);
    Vector<3> v01 = getPixel(x0, y1), v11 = getPixel(x1, y1);
    Vector<3> top = (1.0 - dx) * v00 + dx * v10;
    Vector<3> bottom = (1.0 - dx) * v01 + dx * v11;
    Vector<3> val = ((1.0 - dy) * top + dy * bottom) * sourceIntensity_;

    ScalarType value = extractChannel<ScalarType>(val, channel_);
    if (is_zero)
        *is_zero = scalarL1Norm(value) < EPSILON;
    return value;
}

// ============================================================================
// sampleBoundaryTexture (3D)
// ============================================================================

template<typename ScalarType, int DIM>
ScalarType
ObjScene<ScalarType, DIM>::sampleBoundaryTexture(unsigned char* data, int w, int h, int c, const Vector<2>& uv) const
{
    if (!data || w <= 0 || h <= 0)
        return extractChannel<ScalarType>(Vector<3>(1, 0, 1), channel_);

    double u = uv[0] - std::floor(uv[0]);
    double v = 1.0 - (uv[1] - std::floor(uv[1]));

    double x = u * (w - 1), y = v * (h - 1);
    int x0 = std::max(0, std::min(static_cast<int>(std::floor(x)), w - 1));
    int y0 = std::max(0, std::min(static_cast<int>(std::floor(y)), h - 1));
    int x1 = std::min(x0 + 1, w - 1), y1 = std::min(y0 + 1, h - 1);
    double dx = x - x0, dy = y - y0;

    auto getPixel = [&](int px, int py) -> Vector<3> {
        int idx = (py * w + px) * c;
        double r = data[idx] / 255.0;
        double g = (c > 1) ? data[idx + 1] / 255.0 : r;
        double b = (c > 2) ? data[idx + 2] / 255.0 : 0.0;
        return Vector<3>(r, g, b);
    };

    Vector<3> v00 = getPixel(x0, y0), v10 = getPixel(x1, y0);
    Vector<3> v01 = getPixel(x0, y1), v11 = getPixel(x1, y1);
    Vector<3> top = (1.0 - dx) * v00 + dx * v10;
    Vector<3> bottom = (1.0 - dx) * v01 + dx * v11;
    return extractChannel<ScalarType>((1.0 - dy) * top + dy * bottom, channel_);
}

// ============================================================================
// boundaryCondition
// ============================================================================

template<typename ScalarType, int DIM>
BoundaryCondition<ScalarType>
ObjScene<ScalarType, DIM>::boundaryCondition(const fcpw::Interaction<DIM>& interaction,
                                             const VectorType& queryPoint,
                                             BoundaryType type) const
{
    if constexpr (DIM == 2) {
        auto valFromData = [&](const EdgeBoundaryData& data) -> ScalarType {
            float t = interaction.uv[0];
            ScalarType val = scalarLerp(data.startVal, data.endVal, t);
            if (data.isDoubleSided) {
                auto dir = interaction.p - queryPoint.template cast<float>();
                if (dir.dot(interaction.n) < 0)
                    val = scalarLerp(data.startVal2, data.endVal2, t);
            }
            return val;
        };

        auto normalPointsTowardQuery = [&]() -> bool {
            return (queryPoint.template cast<float>() - interaction.p).dot(interaction.n) > 0;
        };

        if (type == BoundaryType::Dirichlet) {
            return BoundaryCondition<ScalarType>::Dirichlet(valFromData(dirichletData_[interaction.primitiveIndex]));
        } else if (type == BoundaryType::Neumann) {
            ScalarType val = valFromData(neumannData_[interaction.primitiveIndex]);
            if (normalPointsTowardQuery())
                val = scalarNegate(val);
            return BoundaryCondition<ScalarType>::Neumann(val);
        } else if (type == BoundaryType::Robin) {
            ScalarType aVal = valFromData(dirichletData_[interaction.primitiveIndex]);
            ScalarType bVal = valFromData(neumannData_[interaction.primitiveIndex]);
            if (normalPointsTowardQuery())
                bVal = scalarNegate(bVal);
            return BoundaryCondition<ScalarType>(BoundaryType::Robin, aVal, bVal);
        }
    } else {
        int primIdx = interaction.primitiveIndex;
        float w0 = interaction.uv[0], w1 = interaction.uv[1], w2 = 1.0f - w0 - w1;

        auto sampleUV = [&](const std::vector<std::array<Vector<2>, 3>>& uvs,
                            unsigned char* tex,
                            int tw,
                            int th,
                            int tc) -> ScalarType {
            if (primIdx < 0 || primIdx >= static_cast<int>(uvs.size()))
                return ScalarType(0.0);
            const auto& triUVs = uvs[primIdx];
            Vector<2> uv = w0 * triUVs[0] + w1 * triUVs[1] + w2 * triUVs[2];
            return sampleBoundaryTexture(tex, tw, th, tc, uv);
        };

        if (type == BoundaryType::Dirichlet) {
            return BoundaryCondition<ScalarType>::Dirichlet(
              sampleUV(dirichletUVs_, dTexData_, dTexWidth_, dTexHeight_, dTexChannels_));
        } else if (type == BoundaryType::Neumann) {
            ScalarType bVal = sampleUV(neumannUVs_, nTexData_, nTexWidth_, nTexHeight_, nTexChannels_);
            if ((queryPoint.template cast<float>() - interaction.p).dot(interaction.n) > 0)
                bVal = scalarNegate(bVal);
            return BoundaryCondition<ScalarType>::Neumann(bVal);
        } else if (type == BoundaryType::Robin) {
            ScalarType aVal = sampleUV(dirichletUVs_, dTexData_, dTexWidth_, dTexHeight_, dTexChannels_);
            ScalarType bVal = sampleUV(neumannUVs_, nTexData_, nTexWidth_, nTexHeight_, nTexChannels_);
            return BoundaryCondition<ScalarType>(BoundaryType::Robin, aVal, bVal);
        }
    }

    spdlog::error("ObjScene<{}>: Unknown boundary type.", DIM);
    throw std::runtime_error("ObjScene: Unknown boundary type.");
}

// ============================================================================
// generateTargetPoints
// ============================================================================

template<typename ScalarType, int DIM>
std::pair<std::vector<Vector<DIM>>, std::vector<std::size_t>>
ObjScene<ScalarType, DIM>::generateTargetPoints(int width, int height) const
{
    std::vector<Vector<DIM>> points;
    std::vector<std::size_t> pixelIndices;
    points.reserve(static_cast<std::size_t>(width) * height);
    pixelIndices.reserve(static_cast<std::size_t>(width) * height);

    if constexpr (DIM == 2) {
        double dx = (this->renderMax_[0] - this->renderMin_[0]) / width;
        double dy = (this->renderMax_[1] - this->renderMin_[1]) / height;

        for (int y = height - 1; y >= 0; --y) {
            for (int x = 0; x < width; ++x) {
                Vector<DIM> p;
                p[0] = this->renderMin_[0] + (x + 0.5) * dx;
                p[1] = this->renderMin_[1] + (y + 0.5) * dy;
                if (!this->isInsideROI(p))
                    continue;
                points.push_back(p);
                pixelIndices.push_back(static_cast<std::size_t>(height - 1 - y) * width + static_cast<std::size_t>(x));
            }
        }
    } else {
        double dx = solverPlaneWidth_ / width;
        double dy = solverPlaneHeight_ / height;

        for (int y = height - 1; y >= 0; --y) {
            for (int x = 0; x < width; ++x) {
                double u = (x + 0.5) * dx - 0.5 * solverPlaneWidth_;
                double v = (y + 0.5) * dy - 0.5 * solverPlaneHeight_;
                Vector<DIM> p = solverPlaneOrigin_ + u * solverPlaneUAxis_ + v * solverPlaneVAxis_;
                if (!this->isInsideROI(p))
                    continue;
                points.push_back(p);
                pixelIndices.push_back(static_cast<std::size_t>(height - 1 - y) * width + static_cast<std::size_t>(x));
            }
        }
    }

    return { std::move(points), std::move(pixelIndices) };
}

// ============================================================================
// Explicit instantiations
// ============================================================================

template class WOS_SCENE_API ObjScene<Scalar<1>, 2>;
template class WOS_SCENE_API ObjScene<Scalar<3>, 2>;
template class WOS_SCENE_API ObjScene<Scalar<1>, 3>;
template class WOS_SCENE_API ObjScene<Scalar<3>, 3>;

WOS_NAMESPACE_CLOSE_SCOPE
