#pragma once

#include <Eigen/Geometry>

#include "SIMPLib/Geometry/TriangleGeom.h"

#include "DREAM3DReviewFilters/util/TriMesh.h"

/**
 * @brief The Delaunay2D class
 */
class Delaunay2D : public Observable
{
public:
  using Self = Delaunay2D;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  ~Delaunay2D() override;

  static Pointer New(TriMesh::VertexCoordList vertices, double offset, double tolerance, double alpha, Observable* observable = nullptr)
  {
    Pointer sharedPtr(new Delaunay2D(vertices, offset, tolerance, alpha, observable));
    return sharedPtr;
  }

  /**
   * @brief Setter property for Vertices
   */
  void setVertices(const TriMesh::VertexCoordList& value);

  /**
   * @brief Getter property for Vertices
   * @return Value of Vertices
   */
  TriMesh::VertexCoordList getVertices() const;

  /**
   * @brief Setter property for Offset
   */
  void setOffset(const double& value);

  /**
   * @brief Getter property for Offset
   * @return Value of Offset
   */
  double getOffset() const;

  /**
   * @brief Setter property for Tolerance
   */
  void setTolerance(const double& value);

  /**
   * @brief Getter property for Tolerance
   * @return Value of Tolerance
   */
  double getTolerance() const;

  /**
   * @brief Setter property for Alpha
   */
  void setAlpha(const double& value);

  /**
   * @brief Getter property for Alpha
   * @return Value of Alpha
   */
  double getAlpha() const;

  /**
   * @brief Setter property for MessagePrefix
   */
  void setMessagePrefix(const QString& value);

  /**
   * @brief Getter property for MessagePrefix
   * @return Value of MessagePrefix
   */
  QString getMessagePrefix() const;

  /**
   * @brief Setter property for MessageTitle
   */
  void setMessageTitle(const QString& value);

  /**
   * @brief Getter property for MessageTitle
   * @return Value of MessageTitle
   */
  QString getMessageTitle() const;

  /**
   * @brief Setter property for Observer
   */
  void setObserver(Observable* value);

  /**
   * @brief Getter property for Observer
   * @return Value of Observer
   */
  Observable* getObserver() const;

  TriangleGeom::Pointer triangulate();

protected:
  explicit Delaunay2D(TriMesh::VertexCoordList vertices, double offset, double tolerance, double alpha, Observable* observable = nullptr);

private:
  TriMesh::VertexCoordList m_Vertices = {};
  double m_Offset = {};
  double m_Tolerance = {};
  double m_Alpha = {};
  QString m_MessagePrefix = {};
  QString m_MessageTitle = {};
  Observable* m_Observer = nullptr;

  double m_PointBounds[6];
  TriMesh::Pointer m_Delaunay;
  size_t m_NumDuplicatePoints;
  size_t m_NumDegeneracies;

  void initialize();

  int64_t findTriangle(double x[3], int64_t tri, double tol, int64_t nei[3], int64_t pts[3]);

  void checkEdge(int64_t point, double x[3], int64_t p1, int64_t p2, int64_t tri, bool recursive);

  void findPointBounds(TriMesh::VertexCoordList& vertices);

  void fixupBoundaryTriangles(int64_t numVerts, TriMesh::VertexCoordList& points, std::vector<int64_t>& triUse);

  double circumcircle(double a[2], double b[2], double c[2], double center[2]);

  int32_t inCircumcircle(double p[3], double a[3], double b[3], double c[3]);

  double Determinant3x3(double g[3][3]);

  double Dot2D(double a[2], double b[2]);

  double Normalize2x1(double n[2]);

  void findTriangleNormal(int64_t triPoints[3], double n[3]);

  Eigen::Transform<double, 3, Eigen::Affine> findProjectionPlane();

  Delaunay2D(const Delaunay2D&) = delete;     // Copy Constructor Not Implemented
  void operator=(const Delaunay2D&) = delete; // Operator '=' Not Implemented
};
