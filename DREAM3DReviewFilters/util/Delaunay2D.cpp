#include "Delaunay2D.h"

#include <Eigen/Dense>

#include "SIMPLib/Math/MatrixMath.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
Delaunay2D::Delaunay2D(TriMesh::VertexCoordList vertices, double offset, double tolerance, double alpha, Observable* observable)
: m_Vertices(vertices)
, m_Offset(offset)
, m_Tolerance(tolerance)
, m_Alpha(alpha)
, m_Observer(observable)
{
  initialize();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
Delaunay2D::~Delaunay2D() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void Delaunay2D::initialize()
{
  m_PointBounds[0] = std::numeric_limits<double>::max();
  m_PointBounds[1] = std::numeric_limits<double>::lowest();
  m_PointBounds[2] = std::numeric_limits<double>::max();
  m_PointBounds[3] = std::numeric_limits<double>::lowest();
  m_PointBounds[4] = std::numeric_limits<double>::max();
  m_PointBounds[5] = std::numeric_limits<double>::lowest();

  m_NumDuplicatePoints = 0;
  m_NumDegeneracies = 0;

  m_Delaunay = TriMesh::NullPointer();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
TriangleGeom::Pointer Delaunay2D::triangulate()
{
  initialize();

  if(m_Observer != nullptr)
  {
    connect(this, SIGNAL(filterGeneratedMessage(const PipelineMessage&)), m_Observer, SLOT(broadcastPipelineMessage(const PipelineMessage&)));
  }

  /* Eigen::Transform<double, 3, Eigen::Affine> transform = */ findProjectionPlane();

  TriMesh::VertexCoordList projectedVertices(m_Vertices);

  auto numVerts = m_Vertices.size();

  for(size_t v = 0; v < numVerts; v++)
  {
    Eigen::Vector3d point(m_Vertices[v][0], m_Vertices[v][1], m_Vertices[v][2]);
    // Eigen::Vector3d transformedPoint = transform * point;
    // projectedVertices[v][0] = transformedPoint(0);
    // projectedVertices[v][1] = transformedPoint(1);
    // projectedVertices[v][2] = transformedPoint(2);
    projectedVertices[v][0] = point(0);
    projectedVertices[v][1] = point(1);
    projectedVertices[v][2] = point(2);
  }

  findPointBounds(projectedVertices);

  double center[3];
  center[0] = (m_PointBounds[0] + m_PointBounds[1]) / 2.0;
  center[1] = (m_PointBounds[2] + m_PointBounds[3]) / 2.0;
  center[2] = (m_PointBounds[4] + m_PointBounds[5]) / 2.0;

  double diff = 0.0;
  double l = 0.0;

  for(size_t i = 0; i < 3; i++)
  {
    diff = static_cast<double>(m_PointBounds[2 * i + 1]) - static_cast<double>(m_PointBounds[2 * i]);
    l += diff * diff;
  }
  double tol = sqrt(l);
  double radius = m_Offset * tol;
  tol *= m_Tolerance;

  double x[3];

  for(size_t ptId = 0; ptId < 8; ptId++)
  {
    x[0] = center[0] + radius * cos(ptId * 45.0 * SIMPLib::Constants::k_PiOver180D);
    x[1] = center[1] + radius * sin(ptId * 45.0 * SIMPLib::Constants::k_PiOver180D);
    x[2] = center[2];
    projectedVertices.emplace_back(std::initializer_list<float>{float(x[0]), float(x[1]), float(x[2])});
  }

  m_Delaunay = TriMesh::New(projectedVertices);
  m_Delaunay->addTriangle(numVerts + 0, numVerts + 1, numVerts + 2);
  m_Delaunay->addTriangle(numVerts + 2, numVerts + 3, numVerts + 4);
  m_Delaunay->addTriangle(numVerts + 4, numVerts + 5, numVerts + 6);
  m_Delaunay->addTriangle(numVerts + 6, numVerts + 7, numVerts + 0);
  m_Delaunay->addTriangle(numVerts + 0, numVerts + 2, numVerts + 6);
  m_Delaunay->addTriangle(numVerts + 2, numVerts + 4, numVerts + 6);

  int64_t nei[3];
  int64_t neiPts[3];
  int64_t tri[4];
  int64_t pts[3];
  int64_t nodes[4][3];
  int64_t p1;
  int64_t p2;

  tri[0] = 0;

  int64_t progIncrement = static_cast<int64_t>(numVerts / 100);
  int64_t prog = 1;
  int64_t progressInt = 0;
  int64_t counter = 0;

  for(int64_t ptId = 0; ptId < int64_t(numVerts); ptId++)
  {
    x[0] = projectedVertices[ptId][0];
    x[1] = projectedVertices[ptId][1];
    x[2] = projectedVertices[ptId][2];

    nei[0] = (-1); // where we are coming from...nowhere initially

    if((tri[0] = findTriangle(x, tri[0], tol, nei, pts)) >= 0)
    {
      if(nei[0] < 0) // in triangle
      {
        // delete this triangle; create three new triangles
        // first triangle is replaced with one of the new ones

        nodes[0][0] = ptId;
        nodes[0][1] = pts[0];
        nodes[0][2] = pts[1];
        m_Delaunay->removeLinkFromTriangle(pts[2], tri[0]);
        m_Delaunay->replaceTriangleVertices(nodes[0][0], nodes[0][1], nodes[0][2], tri[0]);
        m_Delaunay->addLinkToTriangle(ptId, tri[0]);

        nodes[1][0] = ptId;
        nodes[1][1] = pts[1];
        nodes[1][2] = pts[2];
        tri[1] = m_Delaunay->addTriangle(nodes[1][0], nodes[1][1], nodes[1][2]);

        nodes[2][0] = ptId;
        nodes[2][1] = pts[2];
        nodes[2][2] = pts[0];
        tri[2] = m_Delaunay->addTriangle(nodes[2][0], nodes[2][1], nodes[2][2]);

        // Check edge neighbors for Delaunay criterion. If not satisfied, flip
        // edge diagonal. (This is done recursively.)
        checkEdge(ptId, x, pts[0], pts[1], tri[0], true);
        checkEdge(ptId, x, pts[1], pts[2], tri[1], true);
        checkEdge(ptId, x, pts[2], pts[0], tri[2], true);
      }

      else // on triangle edge
      {
        // update cell list
        m_Delaunay->getTriangleVertices(nei[0], neiPts);
        for(size_t i = 0; i < 3; i++)
        {
          if(neiPts[i] != nei[1] && neiPts[i] != nei[2])
          {
            p1 = neiPts[i];
          }
          if(pts[i] != nei[1] && pts[i] != nei[2])
          {
            p2 = pts[i];
          }
        }

        m_Delaunay->removeLinkFromTriangle(nei[2], tri[0]);
        m_Delaunay->removeLinkFromTriangle(nei[2], nei[0]);

        nodes[0][0] = ptId;
        nodes[0][1] = p2;
        nodes[0][2] = nei[1];
        m_Delaunay->replaceTriangleVertices(nodes[0][0], nodes[0][1], nodes[0][2], tri[0]);

        nodes[1][0] = ptId;
        nodes[1][1] = p1;
        nodes[1][2] = nei[1];
        m_Delaunay->replaceTriangleVertices(nodes[1][0], nodes[1][1], nodes[1][2], nei[0]);

        m_Delaunay->addLinkToTriangle(ptId, tri[0]);
        m_Delaunay->addLinkToTriangle(ptId, nei[0]);

        tri[1] = nei[0];

        nodes[2][0] = ptId;
        nodes[2][1] = p2;
        nodes[2][2] = nei[2];
        tri[2] = m_Delaunay->addTriangle(nodes[2][0], nodes[2][1], nodes[2][2]);

        nodes[3][0] = ptId;
        nodes[3][1] = p1;
        nodes[3][2] = nei[2];
        tri[3] = m_Delaunay->addTriangle(nodes[3][0], nodes[3][1], nodes[3][2]);

        // Check edge neighbors for Delaunay criterion.
        for(size_t i = 0; i < 4; i++)
        {
          checkEdge(ptId, x, nodes[i][1], nodes[i][2], tri[i], true);
        }
      }
    } // if triangle found

    else
    {
      tri[0] = 0; // no triangle found
    }

    if(counter > prog)
    {
      progressInt = static_cast<int64_t>((static_cast<float>(counter) / numVerts) * 100.0f);
      QString ss = m_MessageTitle + QObject::tr(" || %1% Complete").arg(progressInt);
      notifyStatusMessage(ss);
      prog = prog + progIncrement;
    }
    counter++;
  } // for all points

  std::vector<int64_t> triUse(m_Delaunay->getNumberOfTriangles(), 1);

  for(int64_t ptId = numVerts; ptId < int64_t((numVerts + 8)); ptId++)
  {
    std::vector<int64_t> neighborTris = m_Delaunay->getTrianglesToVertex(ptId);
    for(auto&& neighbor : neighborTris)
    {
      triUse[neighbor] = 0;
    }
  }

  ////alpha begin
  // if(m_Alpha > 0.0)
  //{
  //  double alpha2 = m_Alpha * m_Alpha;
  //  double x1[3], x2[3], x3[3];
  //  double xx1[3], xx2[3], xx3[3];
  //  int64_t cellId, numNei, ap1, ap2, neighbor;

  //  TriMesh::Pointer alphaVert = TriMesh::New(projectedVertices);

  //  vtkCellArray *alphaVerts = vtkCellArray::New();
  //  alphaVerts->Allocate(numPoints);
  //  vtkCellArray *alphaLines = vtkCellArray::New();
  //  alphaLines->Allocate(numPoints);

  //  char *pointUse = new char[numPoints + 8];
  //  for(ptId = 0; ptId < (numPoints + 8); ptId++)
  //  {
  //    pointUse[ptId] = 0;
  //  }

  //  //traverse all triangles; evaluating Delaunay criterion
  //  for(i = 0; i < numTriangles; i++)
  //  {
  //    if(triUse[i] == 1)
  //    {
  //      this->Mesh->GetCellPoints(i, npts, triPts);

  //      // if any point is one of the bounding points that was added
  //      // at the beginning of the algorithm, then grab the points
  //      // from the variable "points" (this list has the boundary
  //      // points and the original points have been transformed by the
  //      // input transform).  if none of the points are bounding points,
  //      // then grab the points from the variable "inPoints" so the alpha
  //      // criterion is applied in the nontransformed space.
  //      if(triPts[0]<numPoints && triPts[1]<numPoints && triPts[2]<numPoints)
  //      {
  //        inPoints->GetPoint(triPts[0], x1);
  //        inPoints->GetPoint(triPts[1], x2);
  //        inPoints->GetPoint(triPts[2], x3);
  //      }
  //      else
  //      {
  //        points->GetPoint(triPts[0], x1);
  //        points->GetPoint(triPts[1], x2);
  //        points->GetPoint(triPts[2], x3);
  //      }

  //      // evaluate the alpha criterion in 3D
  //      vtkTriangle::ProjectTo2D(x1, x2, x3, xx1, xx2, xx3);
  //      if(vtkTriangle::Circumcircle(xx1, xx2, xx3, center) > alpha2)
  //      {
  //        triUse[i] = 0;
  //      }
  //      else
  //      {
  //        for(int j = 0; j<3; j++)
  //        {
  //          pointUse[triPts[j]] = 1;
  //        }
  //      }
  //    }//if non-deleted triangle
  //  }//for all triangles

  //  //traverse all edges see whether we need to create some
  //  for(cellId = 0, triangles->InitTraversal();
  //      triangles->GetNextCell(npts, triPts); cellId++)
  //  {
  //    if(!triUse[cellId])
  //    {
  //      for(i = 0; i < npts; i++)
  //      {
  //        ap1 = triPts[i];
  //        ap2 = triPts[(i + 1) % npts];

  //        if(this->BoundingTriangulation || (ap1<numPoints && ap2<numPoints))
  //        {
  //          this->Mesh->GetCellEdgeNeighbors(cellId, ap1, ap2, neighbors);
  //          numNei = neighbors->GetNumberOfIds();

  //          if(numNei < 1 || ((neighbor = neighbors->GetId(0)) > cellId
  //            && !triUse[neighbor]))
  //          {//see whether edge is shorter than Alpha

  //            // same argument as above, if one is a boundary point, get
  //            // it using this->GetPoint() which are transformed points. if
  //            // neither of the points are boundary points, get the from
  //            // inPoints (untransformed points) so alpha comparison done
  //            // untransformed space
  //            if(ap1 < numPoints && ap2 < numPoints)
  //            {
  //              inPoints->GetPoint(ap1, x1);
  //              inPoints->GetPoint(ap2, x2);
  //            }
  //            else
  //            {
  //              this->GetPoint(ap1, x1);
  //              this->GetPoint(ap2, x2);
  //            }
  //            if((vtkMath::Distance2BetweenPoints(x1, x2)*0.25) <= alpha2)
  //            {
  //              pointUse[ap1] = 1; pointUse[ap2] = 1;
  //              pts[0] = ap1;
  //              pts[1] = ap2;
  //              alphaLines->InsertNextCell(2, pts);
  //            }//if passed test
  //          }//test edge
  //        }//if valid edge
  //      }//for all edges of this triangle
  //    }//if triangle not output
  //  }//for all triangles

  //  //traverse all points, create vertices if none used
  //  for(ptId = 0; ptId<(numPoints + 8); ptId++)
  //  {
  //    if((ptId < numPoints || this->BoundingTriangulation)
  //       && !pointUse[ptId])
  //    {
  //      pts[0] = ptId;
  //      alphaVerts->InsertNextCell(1, pts);
  //    }
  //  }

  //  // update output
  //  delete[] pointUse;
  //  output->SetVerts(alphaVerts);
  //  alphaVerts->Delete();
  //  output->SetLines(alphaLines);
  //  alphaLines->Delete();
  //}
  ////alpha end

  fixupBoundaryTriangles(numVerts, projectedVertices, triUse);

  // Remove the 8 bounding triangle vertices
  TriMesh::VertexList delVerts = m_Delaunay->getVertices();
  for(size_t i = numVerts; i < 8; i++)
  {
    delVerts.pop_back();
  }

  // std::vector<int64_t> vertsToRemove(8);
  // std::iota(vertsToRemove.begin(), vertsToRemove.end(), numVerts);
  // m_Delaunay->removeVertices(vertsToRemove);

  size_t numGoodTris = 0;

  for(size_t i = 0; i < triUse.size(); i++)
  {
    if(triUse[i] != 0)
    {
      numGoodTris++;
    }
  }

  std::vector<std::vector<float>> untransformedVerts;
  for(auto i = 0; i < numVerts; i++)
  {
    untransformedVerts.push_back(delVerts[i].vert);
  }

  SharedVertexList::Pointer vertices = TriangleGeom::CreateSharedVertexList(untransformedVerts.size());
  TriangleGeom::Pointer triangles = TriangleGeom::CreateGeometry(numGoodTris, vertices, SIMPL::Geometry::TriangleGeometry);
  float* vertPtr = triangles->getVertexPointer(0);
  size_t* triPtr = triangles->getTriPointer(0);

  TriMesh::TriList triList = m_Delaunay->getTriangles();

  size_t triIter = 0;
  int64_t triVerts[3];

  for(auto i = 0; i < triList.size(); i++)
  {
    if(triUse[i] != 0)
    {
      m_Delaunay->getTriangleVertices(i, triVerts);
      triPtr[3 * triIter + 0] = triVerts[0];
      triPtr[3 * triIter + 1] = triVerts[1];
      triPtr[3 * triIter + 2] = triVerts[2];
      triIter++;
    }
  }

  for(auto i = 0; i < untransformedVerts.size(); i++)
  {
    Eigen::Vector3d point(untransformedVerts[i][0], untransformedVerts[i][1], untransformedVerts[i][2]);
    // Eigen::Vector3d transformedPoint = transform.inverse() * point;
    // vertPtr[3 * i + 0] = transformedPoint(0);
    // vertPtr[3 * i + 1] = transformedPoint(1);
    // vertPtr[3 * i + 2] = transformedPoint(2);
    vertPtr[3 * i + 0] = point(0);
    vertPtr[3 * i + 1] = point(1);
    vertPtr[3 * i + 2] = point(2);
  }

  return triangles;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int64_t Delaunay2D::findTriangle(double x[3], int64_t tri, double tol, int64_t nei[3], int64_t pts[3])
{
  int i, j, ir, ic, inside, i2, i3;
  int64_t newNei;
  double p[3][3], n[2], vp[2], vx[2], dp, minProj;

  m_Delaunay->getTriangleVertices(tri, pts);
  m_Delaunay->getVertexCoordinates(pts[0], p[0]);
  m_Delaunay->getVertexCoordinates(pts[1], p[1]);
  m_Delaunay->getVertexCoordinates(pts[2], p[2]);

  srand(tri);
  ir = rand() % 3;
  const double del2D_tolerance = 1.0e-014;

  for(inside = 1, minProj = del2D_tolerance, ic = 0; ic < 3; ic++)
  {
    i = (ir + ic) % 3;
    i2 = (i + 1) % 3;
    i3 = (i + 2) % 3;

    // create a 2D edge normal to define a "half-space"; evaluate points (i.e.,
    // candidate point and other triangle vertex not on this edge).
    n[0] = -(p[i2][1] - p[i][1]);
    n[1] = p[i2][0] - p[i][0];
    Normalize2x1(n);

    // compute local vectors
    for(j = 0; j < 2; j++)
    {
      vp[j] = p[i3][j] - p[i][j];
      vx[j] = x[j] - p[i][j];
    }

    // check for duplicate point
    Normalize2x1(vp);
    if(Normalize2x1(vx) <= tol)
    {
      m_NumDuplicatePoints++;
      return -1;
    }

    // see if two points are in opposite half spaces
    dp = Dot2D(n, vx) * (Dot2D(n, vp) < 0 ? -1.0 : 1.0);
    if(dp < del2D_tolerance)
    {
      if(dp < minProj) // track edge most orthogonal to point direction
      {
        inside = 0;
        nei[1] = pts[i];
        nei[2] = pts[i2];
        minProj = dp;
      }
    } // outside this edge
  }   // for each edge

  if(inside != 0) // all edges have tested positive
  {
    nei[0] = (-1);
    return tri;
  }

  if(!inside && (fabs(minProj) < del2D_tolerance)) // on edge
  {
    nei[0] = m_Delaunay->getTriangleEdgeNeighbor(nei[1], nei[2], tri);
    return tri;
  }

  else // walk towards point
  {
    newNei = m_Delaunay->getTriangleEdgeNeighbor(nei[1], nei[2], tri);
    if(newNei == nei[0])
    {
      m_NumDegeneracies++;
      return -1;
    }
    else
    {
      nei[0] = tri;
      return findTriangle(x, newNei, tol, nei, pts);
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void Delaunay2D::checkEdge(int64_t point, double x[3], int64_t p1, int64_t p2, int64_t tri, bool recursive)
{
  double x1[3], x2[3], x3[3];
  int64_t neighbor;

  m_Delaunay->getVertexCoordinates(p1, x1);
  m_Delaunay->getVertexCoordinates(p2, x2);

  neighbor = m_Delaunay->getTriangleEdgeNeighbor(p1, p2, tri);

  if(neighbor > 0) // i.e., not a boundary edge
  {
    // get neighbor info including opposite point
    int64_t oppositeVert = m_Delaunay->getOppositeVertex(p1, p2, neighbor);
    if(oppositeVert < 0)
    {
      return;
    }
    m_Delaunay->getVertexCoordinates(oppositeVert, x3);

    // see whether point is in circumcircle
    if(inCircumcircle(x3, x, x1, x2) != 0)
    { // swap diagonal
      m_Delaunay->removeLinkFromTriangle(p1, tri);
      m_Delaunay->removeLinkFromTriangle(p2, neighbor);
      m_Delaunay->addLinkToTriangle(point, neighbor);
      m_Delaunay->addLinkToTriangle(oppositeVert, tri);

      m_Delaunay->replaceTriangleVertices(point, oppositeVert, p2, tri);
      m_Delaunay->replaceTriangleVertices(point, p1, oppositeVert, neighbor);

      if(recursive)
      {
        // two new edges become suspect
        checkEdge(point, x, oppositeVert, p2, tri, true);
        checkEdge(point, x, p1, oppositeVert, neighbor, true);
      }
    } // in circle
  }   // interior edge
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void Delaunay2D::findPointBounds(TriMesh::VertexCoordList& vertices)
{
  for(auto&& vert : vertices)
  {
    if(vert[0] < m_PointBounds[0])
    {
      m_PointBounds[0] = static_cast<double>(vert[0]);
    }
    if(vert[0] > m_PointBounds[1])
    {
      m_PointBounds[1] = static_cast<double>(vert[0]);
    }
    if(vert[1] < m_PointBounds[2])
    {
      m_PointBounds[2] = static_cast<double>(vert[1]);
    }
    if(vert[1] > m_PointBounds[3])
    {
      m_PointBounds[3] = static_cast<double>(vert[1]);
    }
    if(vert[2] < m_PointBounds[4])
    {
      m_PointBounds[4] = static_cast<double>(vert[2]);
    }
    if(vert[2] > m_PointBounds[5])
    {
      m_PointBounds[5] = static_cast<double>(vert[2]);
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void Delaunay2D::fixupBoundaryTriangles(int64_t numVerts, TriMesh::VertexCoordList& points, std::vector<int64_t>& triUse)
{
  bool isConnected;
  int64_t numSwaps = 0;
  int64_t p1;
  int64_t p2;
  int64_t p3;
  int64_t triPts[3];
  int64_t neiPts[3];
  int64_t pts[3];
  int64_t swapPts[3];
  double n1[3];
  double n2[3];

  for(int64_t ptId = 0; ptId < numVerts; ptId++)
  {
    // check if point is only connected to triangles scheduled for
    // removal
    std::vector<int64_t> neighbors = m_Delaunay->getTrianglesToVertex(ptId);
    auto ncells = neighbors.size();

    isConnected = false;

    for(size_t i = 0; i < ncells; i++)
    {
      if(triUse[neighbors[i]] != 0)
      {
        isConnected = true;
        break;
      }
    }

    // this point will be connected in the output
    if(isConnected)
    {
      // point is connected: continue
      continue;
    }

    // This point is only connected to triangles scheduled for removal.
    // Therefore it will not be connected in the output triangulation.
    // Let's swap edges to create a triangle with 3 inner points.
    // - inner points have an id < numPoints
    // - boundary point ids are, numPoints <= id < numPoints+8.

    // visit every edge connected to that point.
    // check the 2 triangles touching at that edge.
    // if one triangle is connected to 2 non-boundary points

    for(auto i = 0; i < ncells; i++)
    {
      int64_t tri1 = neighbors[i];
      m_Delaunay->getTriangleVertices(tri1, triPts);

      if(triPts[0] == ptId)
      {
        p1 = triPts[1];
        p2 = triPts[2];
      }
      else if(triPts[1] == ptId)
      {
        p1 = triPts[2];
        p2 = triPts[0];
      }
      else
      {
        p1 = triPts[0];
        p2 = triPts[1];
      }

      // if both p1 & p2 are boundary points,
      // we skip them.
      if(p1 >= numVerts && p2 >= numVerts)
      {
        continue;
      }

      int64_t tri2 = m_Delaunay->getTriangleEdgeNeighbor(p1, p2, tri1);

      // get the 3 points of the neighbor triangle
      m_Delaunay->getTriangleVertices(tri2, neiPts);

      // locate the point different from p1 and p2
      if(neiPts[0] != p1 && neiPts[0] != p2)
      {
        p3 = neiPts[0];
      }
      else if(neiPts[1] != p1 && neiPts[1] != p2)
      {
        p3 = neiPts[1];
      }
      else
      {
        p3 = neiPts[2];
      }

      // create the two new triangles.
      // we just need to replace their pt ids.
      pts[0] = ptId;
      pts[1] = p1;
      pts[2] = p3;

      swapPts[0] = ptId;
      swapPts[1] = p3;
      swapPts[2] = p2;

      findTriangleNormal(pts, n1);
      findTriangleNormal(swapPts, n2);

      // the normals must be along the same direction,
      // or one triangle is upside down.
      if(MatrixMath::DotProduct3x1(n1, n2) < 0.0)
      {
        // do not swap diagonal
        continue;
      }

      // swap edge [p1 p2] and diagonal [ptId p3]

      // it's ok to swap the diagonal
      m_Delaunay->removeLinkFromTriangle(p1, tri2);
      m_Delaunay->removeLinkFromTriangle(p2, tri1);
      m_Delaunay->addLinkToTriangle(ptId, tri2);
      m_Delaunay->addLinkToTriangle(p3, tri1);

      m_Delaunay->replaceTriangleVertices(pts[0], pts[1], pts[2], tri1);
      m_Delaunay->replaceTriangleVertices(swapPts[0], swapPts[1], swapPts[2], tri2);

      triUse[tri1] = (p1 < numVerts && p3 < numVerts);
      triUse[tri2] = (p3 < numVerts && p2 < numVerts);

      // update the 'scheduled for removal' flag of the first triangle.
      // The second triangle was not scheduled for removal anyway.
      numSwaps++;
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
double Delaunay2D::circumcircle(double a[2], double b[2], double c[2], double center[2])
{
  double n12[2] = {0.0, 0.0};
  double n13[2] = {0.0, 0.0};
  double x12[2] = {0.0, 0.0};
  double x13[2] = {0.0, 0.0};
  double sum = 0.0;
  double diff = 0.0;
  size_t i;

  for(i = 0; i < 2; i++)
  {
    n12[i] = b[i] - a[i];
    n13[i] = c[i] - a[i];
    x12[i] = (b[i] + a[i]) / 2.0;
    x13[i] = (c[i] + a[i]) / 2.0;
  }

  Eigen::Matrix2d A;
  A << n12[0], n12[1], n13[0], n13[1];
  Eigen::Vector2d rhs;
  rhs << Dot2D(n12, x12), Dot2D(n13, x13);
  Eigen::Vector2d solution = A.colPivHouseholderQr().solve(rhs);

  center[0] = solution(0);
  center[1] = solution(1);

  for(sum = 0, i = 0; i < 2; i++)
  {
    diff = a[i] - center[i];
    sum += diff * diff;
    diff = b[i] - center[i];
    sum += diff * diff;
    diff = c[i] - center[i];
    sum += diff * diff;
  }

  if((sum /= 3.0) > std::numeric_limits<double>::max())
  {
    return std::numeric_limits<double>::max();
  }

  return sum;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int32_t Delaunay2D::inCircumcircle(double p[3], double a[3], double b[3], double c[3])
{
  double radius2 = 0.0;
  double center[2] = {0.0, 0.0};
  double dist2 = 0.0;

  radius2 = circumcircle(a, b, c, center);

  dist2 = (p[0] - center[0]) * (p[0] - center[0]) + (p[1] - center[1]) * (p[1] - center[1]);

  if(dist2 < (0.999999999999 * radius2))
  {
    return 1;
  }

  return 0;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
double Delaunay2D::Determinant3x3(double g[3][3])
{
  return (g[0][0] * (g[1][1] * g[2][2] - g[1][2] * g[2][1])) - (g[0][1] * (g[1][0] * g[2][2] - g[1][2] * g[2][0])) + (g[0][2] * (g[1][0] * g[2][1] - g[1][1] * g[2][0]));
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
double Delaunay2D::Dot2D(double a[2], double b[2])
{
  return (a[0] * b[0] + a[1] * b[1]);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
double Delaunay2D::Normalize2x1(double n[2])
{
  double denom;
  denom = sqrt(((n[0] * n[0]) + (n[1] * n[1])));
  if(denom != 0)
  {
    n[0] /= denom;
    n[1] /= denom;
  }
  return denom;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void Delaunay2D::findTriangleNormal(int64_t triPoints[3], double n[3])
{
  double length;

  double ax, ay, az, bx, by, bz;

  double a[3];
  double b[3];
  double c[3];

  m_Delaunay->getVertexCoordinates(triPoints[0], a);
  m_Delaunay->getVertexCoordinates(triPoints[1], b);
  m_Delaunay->getVertexCoordinates(triPoints[2], c);

  // order is important!!! maintain consistency with triangle vertex order
  ax = c[0] - b[0];
  ay = c[1] - b[1];
  az = c[2] - b[2];
  bx = a[0] - b[0];
  by = a[1] - b[1];
  bz = a[2] - b[2];

  n[0] = (ay * bz - az * by);
  n[1] = (az * bx - ax * bz);
  n[2] = (ax * by - ay * bx);

  if((length = sqrt((n[0] * n[0] + n[1] * n[1] + n[2] * n[2]))) != 0.0)
  {
    n[0] /= length;
    n[1] /= length;
    n[2] /= length;
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
Eigen::Transform<double, 3, Eigen::Affine> Delaunay2D::findProjectionPlane()
{
  auto numVertices = m_Vertices.size();
  double matrix[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  double normal[3] = {0.0, 0.0, 0.0};
  double origin[3] = {0.0, 0.0, 0.0};
  double v[3] = {0.0, 0.0, 0.0};

  double m[9];
  double *c1, *c2, *c3, det;

  const double tolerance = 1.0e-03;

  for(auto&& vert : m_Vertices)
  {
    v[0] += vert[0] * vert[2];
    v[1] += vert[1] * vert[2];
    v[2] += vert[2];

    matrix[0][0] += vert[0] * vert[0];
    matrix[0][1] += vert[0] * vert[1];
    matrix[0][2] += vert[0];

    matrix[1][0] += vert[0] * vert[1];
    matrix[1][1] += vert[1] * vert[1];
    matrix[1][2] += vert[1];

    matrix[2][0] += vert[0];
    matrix[2][1] += vert[1];
  }

  matrix[2][2] = numVertices;

  origin[0] = matrix[0][2] / numVertices;
  origin[1] = matrix[1][2] / numVertices;
  origin[2] = v[2] / numVertices;

  c1 = m;
  c2 = m + 3;
  c3 = m + 6;

  double tmpMatrixInv[3][3] = {{matrix[0][0], matrix[1][0], matrix[2][0]}, {matrix[0][1], matrix[1][1], matrix[2][1]}, {matrix[0][2], matrix[1][2], matrix[2][2]}};

  double tmpMatrix0[3][3] = {{v[0], matrix[1][0], matrix[2][0]}, {v[1], matrix[1][1], matrix[2][1]}, {v[2], matrix[1][2], matrix[2][2]}};

  double tmpMatrix1[3][3] = {{matrix[0][0], v[0], matrix[2][0]}, {matrix[0][1], v[1], matrix[2][1]}, {matrix[0][2], v[2], matrix[2][2]}};

  if((det = Determinant3x3(tmpMatrixInv)) > tolerance)
  {
    normal[0] = Determinant3x3(tmpMatrix0) / det;
    normal[1] = Determinant3x3(tmpMatrix1) / det;
    normal[2] = -1.0;
  }

  double zAxis[3] = {0.0, 0.0, 1.0};
  double rotationAxis[3] = {0.0, 0.0, 0.0};

  MatrixMath::Normalize3x1(normal);
  MatrixMath::CrossProduct(normal, zAxis, rotationAxis);
  MatrixMath::Normalize3x1(rotationAxis);

  double rotationAngle = acos(MatrixMath::DotProduct3x1(zAxis, normal));

  Eigen::Transform<double, 3, Eigen::Affine> transform =
      Eigen::AngleAxisd(rotationAngle, Eigen::Vector3d(rotationAxis[0], rotationAxis[1], rotationAxis[2])) * Eigen::Translation3d(-origin[0], -origin[1], -origin[2]);

  return transform;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
Delaunay2D::Pointer Delaunay2D::NullPointer()
{
  return Pointer(static_cast<Self*>(nullptr));
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void Delaunay2D::setVertices(const TriMesh::VertexCoordList& value)
{
  m_Vertices = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
TriMesh::VertexCoordList Delaunay2D::getVertices() const
{
  return m_Vertices;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void Delaunay2D::setOffset(const double& value)
{
  m_Offset = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
double Delaunay2D::getOffset() const
{
  return m_Offset;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void Delaunay2D::setTolerance(const double& value)
{
  m_Tolerance = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
double Delaunay2D::getTolerance() const
{
  return m_Tolerance;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void Delaunay2D::setAlpha(const double& value)
{
  m_Alpha = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
double Delaunay2D::getAlpha() const
{
  return m_Alpha;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void Delaunay2D::setMessagePrefix(const QString& value)
{
  m_MessagePrefix = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString Delaunay2D::getMessagePrefix() const
{
  return m_MessagePrefix;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void Delaunay2D::setMessageTitle(const QString& value)
{
  m_MessageTitle = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
QString Delaunay2D::getMessageTitle() const
{
  return m_MessageTitle;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void Delaunay2D::setObserver(Observable* value)
{
  m_Observer = value;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
Observable* Delaunay2D::getObserver() const
{
  return m_Observer;
}
