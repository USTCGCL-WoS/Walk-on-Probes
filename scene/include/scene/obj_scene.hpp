#pragma once

#include <array>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "core/scalar.hpp"
#include "scene/poisson.hpp"
#include "scene/scene_api.h"

WOS_NAMESPACE_OPEN_SCOPE

enum class ChannelMode
{
    Red,
    Green,
    Blue,
    RGB
};

template<typename ScalarType, int DIM>
class ObjScene : public PoissonScene<ScalarType, DIM>
{
  public:
    using VectorType = Vector<DIM>;

    explicit ObjScene(const std::string& sceneDir, ChannelMode channel = ChannelMode::RGB);
    virtual ~ObjScene();

    ScalarType source(const VectorType& p) const override;
    ScalarType source(const VectorType& p, bool* is_zero) const;

    BoundaryCondition<ScalarType> boundaryCondition(const fcpw::Interaction<DIM>& interaction,
                                                    const VectorType& queryPoint,
                                                    BoundaryType type) const override;

    std::pair<std::vector<Vector<DIM>>, std::vector<std::size_t>>
    generateTargetPoints(int width, int height) const override;

  private:
    VectorType sourceMin_;
    VectorType sourceMax_;
    int sourceResolution_;
    ChannelMode channel_;
    std::vector<fcpw::Scene<DIM>*> active_scenes_;

    double sourceIntensity_ = 100.0;
    int texWidth_ = 0, texHeight_ = 0, texChannels_ = 0;
    unsigned char* texData_ = nullptr;

    struct EdgeBoundaryData
    {
        bool isDoubleSided = false;
        ScalarType startVal, endVal, startVal2, endVal2;
    };
    std::vector<EdgeBoundaryData> dirichletData_;
    std::vector<EdgeBoundaryData> neumannData_;

    int dTexWidth_ = 0, dTexHeight_ = 0, dTexChannels_ = 0;
    unsigned char* dTexData_ = nullptr;
    int nTexWidth_ = 0, nTexHeight_ = 0, nTexChannels_ = 0;
    unsigned char* nTexData_ = nullptr;

    std::map<int, int> faceTypes_;
    std::vector<std::array<Vector<2>, 3>> dirichletUVs_;
    std::vector<std::array<Vector<2>, 3>> neumannUVs_;

    void loadConfig(const std::string& path);
    void loadGeometryAndBoundary(const std::string& objPath, const std::string& boundaryPath);
    void loadInterest(const std::string& path);
    void loadSource(const std::string& path);
    void loadBoundaryTextures(const std::string& dir);

    ScalarType sampleTexture(const VectorType& p, bool* is_zero) const;
    ScalarType sampleBoundaryTexture(unsigned char* data, int w, int h, int c, const Vector<2>& uv) const;

    // solver plane (3D)
    Eigen::Vector3d solverPlaneOrigin_{ 0.0, 0.0, 0.0 };
    Eigen::Vector3d solverPlaneUAxis_{ 1.0, 0.0, 0.0 };
    Eigen::Vector3d solverPlaneVAxis_{ 0.0, 0.0, 1.0 };
    double solverPlaneWidth_ = 2.0;
    double solverPlaneHeight_ = 2.0;
};

extern template class ObjScene<Scalar<1>, 2>;
extern template class ObjScene<Scalar<3>, 2>;
extern template class ObjScene<Scalar<1>, 3>;
extern template class ObjScene<Scalar<3>, 3>;

WOS_NAMESPACE_CLOSE_SCOPE
