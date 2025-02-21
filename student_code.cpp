#include "student_code.h"
#include "mutablePriorityQueue.h"

using namespace std;

namespace CGL
{

  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (class member).
   *
   * @param points A vector of points in 2D
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector2D> BezierCurve::evaluateStep(std::vector<Vector2D> const &points)
  { 
    // TODO Part 1.
    // return std::vector<Vector2D>();

    //problem 1 task 1:
      std::vector<Vector2D> nextPoints;
      // Iterate through the control points, stopping one before the end.
      for (size_t i = 0; i < points.size() - 1; i++) {
          // Perform linear interpolation between points[i] and points[i+1].
          Vector2D interpolated = (1 - t) * points[i] + t * points[i + 1];
          nextPoints.push_back(interpolated);
      }
      return nextPoints;
  }
  


  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (function parameter).
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector3D> BezierPatch::evaluateStep(std::vector<Vector3D> const &points, double t) const
  {
    // TODO Part 2.
    
    // return std::vector<Vector3D>();

    // problem 1 task 2
    if (points.size() < 2) {
      std::cerr << "evaluateStep: points.size() < 2, returning points. Size = " << points.size() 
                << ", t = " << t << std::endl;
      return points;
  }
  std::vector<Vector3D> nextPoints;
  for (size_t i = 0; i < points.size() - 1; i++) {
      Vector3D interpolated = (1 - t) * points[i] + t * points[i + 1];
      nextPoints.push_back(interpolated);
  }
  std::cerr << "evaluateStep: reduced from " << points.size() << " to " << nextPoints.size() 
            << ", t = " << t << std::endl;
  return nextPoints;
}


  /**
   * Fully evaluates de Casteljau's algorithm for a vector of points at scalar parameter t
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate1D(std::vector<Vector3D> const &points, double t) const
  {
    // // TODO Part 2.
    // return Vector3D();

    // problem 1 task 2
    std::cerr << "evaluate1D: points.size() = " << points.size() << ", t = " << t << std::endl;
    if (points.empty()) {
        std::cerr << "Error: evaluate1D called with an empty vector of points." << std::endl;
        return Vector3D();
    }
    if (points.size() == 1)
        return points[0];
    std::vector<Vector3D> nextPoints = evaluateStep(points, t);
    return evaluate1D(nextPoints, t);
  }

  /**
   * Evaluates the Bezier patch at parameter (u, v)
   *
   * @param u         Scalar interpolation parameter
   * @param v         Scalar interpolation parameter (along the other axis)
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate(double u, double v) const 
  {  
    // // TODO Part 2.
    // return Vector3D();

    // problem 1 task 2
    std::vector<Vector3D> intermediatePoints;
    for (size_t i = 0; i < controlPoints.size(); i++) {
        if (controlPoints[i].empty()) {
            std::cerr << "Error: controlPoints[" << i << "] is empty." << std::endl;
            continue;
        }
        Vector3D p = evaluate1D(controlPoints[i], u);
        intermediatePoints.push_back(p);
    }
    if (intermediatePoints.empty()) {
        std::cerr << "Error: intermediatePoints is empty." << std::endl;
        return Vector3D();
    }
    return evaluate1D(intermediatePoints, v);
}


  // Vector3D Vertex::normal( void ) const
  // {
    // TODO Part 3.
    // Returns an approximate unit normal at this vertex, computed by
    // taking the area-weighted average of the normals of neighboring
    // triangles, then normalizing.
  //   return Vector3D();
  // }

  // problem 1 task 3
  Vector3D Vertex::normal() const {
    // Accumulate the area-weighted normals.
    Vector3D areaWeightedNormal(0, 0, 0);
    
    // Get the starting halfedge iterator.
    HalfedgeCIter h = this->halfedge();
    // Instead of "if (!h)", we check if h equals a default-constructed iterator.
    if (h == HalfedgeCIter()) return areaWeightedNormal;
    
    // Save the starting halfedge to know when we've looped around.
    HalfedgeCIter start = h;
    
    do {
      // Get the face iterator for the current halfedge.
      FaceCIter f = h->face();
      // Check that the face is valid and not a boundary face.
      if (f != FaceCIter() && !((*f).isBoundary())) {
        // p0 is the position at this vertex.
        Vector3D p0 = this->position;
        // p1 is the vertex from the next halfedge.
        Vector3D p1 = h->next()->vertex()->position;
        // p2 is the vertex from the halfedge after that.
        Vector3D p2 = h->next()->next()->vertex()->position;
        
        // Compute two edge vectors.
        Vector3D e1 = p1 - p0;
        Vector3D e2 = p2 - p0;
        
        // The cross product gives a vector perpendicular to the face.
        Vector3D faceNormal = cross(e1, e2);
        
        // Add the area-weighted normal (note: cross product magnitude is 2x the area).
        areaWeightedNormal += faceNormal;
      }
      // Move to the next halfedge around the vertex.
      h = h->twin()->next();
    } while (h != start);
    
    // Normalize the accumulated vector (using the in-place normalize() method).
    Vector3D normal = areaWeightedNormal;
    normal.normalize();
    return normal;
  }
  

  


  // EdgeIter HalfedgeMesh::flipEdge( EdgeIter e0 )
  // {
  // //   // TODO Part 4.
  // //   // This method should flip the given edge and return an iterator to the flipped edge.
  // //   return EdgeIter();
  // // }

    // problem 1 task 4

    EdgeIter HalfedgeMesh::flipEdge( EdgeIter e0 ) {
      // Do not flip a boundary edge.
      if ( e0->isBoundary() ) return e0;
      
      // Retrieve the two halfedges along the shared edge.
      HalfedgeIter h0  = e0->halfedge();   // h0: from vertex b to vertex c.
      HalfedgeIter h0t = h0->twin();         // h0t: from vertex c to vertex b.
      
      // Do not flip if either adjacent face is a boundary.
      if ( h0->face()->isBoundary() || h0t->face()->isBoundary() )
        return e0;
      
      // below is the current configuration of the triangle and edge
      // the current configuration is a triangle with vertices a, b, c and an edge between b and c
      // the edge is shared by two triangles T₀ = (a, b, c) and T₁ = (c, b, d)

      ////////////////////////////////////
      //       c                c       //  
      //       |                        //
      //  a    |    d -->  a---------d  //  
      //       |                        //
      //       b                b       //  
      ////////////////////////////////////

      // triangle T0 = (a, b, c):
      HalfedgeIter h1 = h0->next();         // h1: from c to a.  => h1->vertex() = a.
      HalfedgeIter h2 = h1->next();         // h2: from a to b.  => h2->vertex() = b.

      // triangle T1 = (c, b, d):
      HalfedgeIter h4 = h0t->next();        // h4: from b to d.  => h4->vertex() = d.
      HalfedgeIter h5 = h4->next();         // h5: from d to c.  => h5->vertex() = c.

      // four key vertices.
      // a = h1->vertex(), b = h2->vertex(), c = h5->vertex(), d = h4->vertex()
      VertexIter c = h1->vertex(); //c,a
      VertexIter a = h2->vertex(); // a,b
      VertexIter d = h5->vertex(); // d,c
      VertexIter b = h4->vertex(); // b,d

      // Flip the shared edge
      // After the flip, the new common edge should connect a and d.
      // Reassign the endpoints of the shared halfedges:
      h0->vertex()  = a;;  // h0 goes from b to d.
      h0t->vertex() = d;;  // h0t goes from c to a.
      // h1->vertex() = d;
      // h5->vertex() = b;


      // Rebuild the face cycles by reassigning the "next" pointers ---
      // For the new triangle T0 (should be (a, d, c)):
      // Use the cycle: (h1, h0, h5).
      h1->next() = h0;
      h0->next() = h5;
      h5->next() = h1;

      // // For the new triangle T1 (should be (a, b, d)):
      // // Use the cycle: (h0t, h2, h4).
      h0t->next() = h2;
      h2->next() = h4;
      h4->next() = h0t;


  
      // // Update face pointers ---
      FaceIter f0 = h0->face();    // Face of T0 remains f0.
      FaceIter f1 = h0t->face();   // Face of T1 remains f1.
      f0->halfedge() = h1;         // Choose h1 as the representative for f0.
      f1->halfedge() = h0t;        // Choose h0t for f1.
      
      h1->face() = f0;  h0->face() = f0;  h5->face() = f0;
      h0t->face() = f1; h2->face() = f1;  h4->face() = f1;


      // Update vertex halfedge pointers 
      a->halfedge() = h1; 
      b->halfedge() = h2; 
      c->halfedge() = h5; 
      d->halfedge() = h4; 

      return e0;


    }

    
    

    


VertexIter HalfedgeMesh::splitEdge( EdgeIter e0 ) {
  // Only process interior edges.
  if ( e0->isBoundary() ) return VertexIter();

  // Retrieve the two halfedges of the shared edge.
  HalfedgeIter h0  = e0->halfedge();   // h0: originally from b to c.
  HalfedgeIter h0t = h0->twin();         // h0t: originally from c to b.
  
  // Reject if either adjacent face is a boundary.
  if ( h0->face()->isBoundary() || h0t->face()->isBoundary() )
    return VertexIter();
  
  // --- Identify original vertices using our convention ---
  // In T₀ = (a, b, c):
  HalfedgeIter h1 = h0->next();         // h1: from c to a → a = h1->vertex()
  HalfedgeIter h2 = h1->next();         // h2: from a to b → b = h2->vertex()
  // In T₁ = (c, b, d):
  HalfedgeIter h4 = h0t->next();        // h4: from b to d → d = h4->vertex()
  HalfedgeIter h5 = h4->next();         // h5: from d to c → c = h5->vertex()
  
  VertexIter a = h1->vertex();
  VertexIter b = h2->vertex();
  VertexIter d = h4->vertex();
  VertexIter c = h5->vertex();
  
  // --- Create new vertex m at the midpoint of edge (b, c) ---
  VertexIter m = newVertex();
  m->position = 0.5 * (b->position + c->position);
  
  // --- SPLIT THE SHARED EDGE ---
  // We split edge (b, c) into two:
  // Reuse e0 for edge (b, m): update h0 so that it now goes from b to m.
  h0->vertex() = m;
  // Reuse the original twin for the twin halfedge: set h0t->vertex() = b, so h0t becomes m → b.
  h0t->vertex() = b;
  // Create a new edge for (m, c).
  EdgeIter e_new = newEdge();
  HalfedgeIter h_new = newHalfedge();    // h_new: from m to c.
  HalfedgeIter h_new_t = newHalfedge();   // h_new_t: from c to m.
  h_new->vertex() = c;
  h_new_t->vertex() = m;
  h_new->twin() = h_new_t;
  h_new_t->twin() = h_new;
  e_new->halfedge() = h_new;
  h_new->edge() = e_new;
  h_new_t->edge() = e_new;
  
  // --- INSERT NEW CONNECTING EDGES ---
  // Create edge from m to a.
  EdgeIter e_AM = newEdge();
  HalfedgeIter h_AM   = newHalfedge();    // from m to a.
  HalfedgeIter h_AM_t = newHalfedge();    // from a to m.
  h_AM->vertex() = a;
  h_AM_t->vertex() = m;
  h_AM->twin() = h_AM_t;
  h_AM_t->twin() = h_AM;
  e_AM->halfedge() = h_AM;
  h_AM->edge() = e_AM;
  h_AM_t->edge() = e_AM;
  
  // Create edge from m to d.
  EdgeIter e_DM = newEdge();
  HalfedgeIter h_DM   = newHalfedge();    // from m to d.
  HalfedgeIter h_DM_t = newHalfedge();    // from d to m.
  h_DM->vertex() = d;
  h_DM_t->vertex() = m;
  h_DM->twin() = h_DM_t;
  h_DM_t->twin() = h_DM;
  e_DM->halfedge() = h_DM;
  h_DM->edge() = e_DM;
  h_DM_t->edge() = e_DM;
  
  // --- REBUILD FACE CYCLES ---
  // We will now subdivide the original two triangles into four new faces.
  
  // Face F0: (a, b, m) from T₀.
  // Cycle: edge from a to b: use h2 (a→b), then edge from b to m: use h0 (b→m), then edge from m to a: use h_AM (m→a).
  FaceIter f0 = h0->face();  // reuse T₀ for F0.
  h2->next() = h0;
  h0->next() = h_AM;
  h_AM->next() = h2;
  f0->halfedge() = h2;
  h2->face() = f0;
  h0->face() = f0;
  h_AM->face() = f0;
  
  // Face F1: (a, m, c) from T₀.
  // Cycle: edge from a to m: use h_AM_t (a→m), then edge from m to c: use h_new (m→c), then edge from c to a: use h1 (c→a).
  FaceIter f1 = newFace();
  h_AM_t->next() = h_new;
  h_new->next() = h1;
  h1->next() = h_AM_t;
  f1->halfedge() = h_AM_t;
  h_AM_t->face() = f1;
  h_new->face() = f1;
  h1->face() = f1;
  
  // Face F2: (d, c, m) from T₁.
  // Cycle: edge from d to c: use h5 (d→c), then edge from c to m: use h_new_t (c→m), then edge from m to d: use h_DM (m→d).
  FaceIter f2 = h0t->face();  // reuse T₁ for F2.
  h5->next() = h_new_t;
  h_new_t->next() = h_DM;
  h_DM->next() = h5;
  f2->halfedge() = h5;
  h5->face() = f2;
  h_new_t->face() = f2;
  h_DM->face() = f2;
  
  // Face F3: (d, m, b) from T₁.
  // Cycle: edge from d to m: use h_DM_t (d→m), then edge from m to b: use h0t (m→b), then edge from b to d: use h4 (b→d).
  FaceIter f3 = newFace();
  h_DM_t->next() = h0t;
  h0t->next() = h4;
  h4->next() = h_DM_t;
  f3->halfedge() = h_DM_t;
  h_DM_t->face() = f3;
  h0t->face() = f3;
  h4->face() = f3;
  
  // --- UPDATE VERTEX HALFEDGE POINTERS (arbitrarily) ---
  a->halfedge() = h1;
  b->halfedge() = h2;
  c->halfedge() = h_new;   // or h_new_t
  d->halfedge() = h4;
  m->halfedge() = h_AM;    // for example.
  
  return m;
}




//   void MeshResampler::upsample( HalfedgeMesh& mesh )
//   {
//     // TODO Part 6.
//     // This routine should increase the number of triangles in the mesh using Loop subdivision.
//     // One possible solution is to break up the method as listed below.

//     // 1. Compute new positions for all the vertices in the input mesh, using the Loop subdivision rule,
//     // and store them in Vertex::newPosition. At this point, we also want to mark each vertex as being
//     // a vertex of the original mesh.
    
//     // 2. Compute the updated vertex positions associated with edges, and store it in Edge::newPosition.
    
//     // 3. Split every edge in the mesh, in any order. For future reference, we're also going to store some
//     // information about which subdivide edges come from splitting an edge in the original mesh, and which edges
//     // are new, by setting the flat Edge::isNew. Note that in this loop, we only want to iterate over edges of
//     // the original mesh---otherwise, we'll end up splitting edges that we just split (and the loop will never end!)
    
//     // 4. Flip any new edge that connects an old and new vertex.

//     // 5. Copy the new vertex positions into final Vertex::position.

//   }
// }

void MeshResampler::upsample( HalfedgeMesh &mesh ) {
  // Step 0: Mark all current vertices and edges as old.
  for (VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++) {
    v->isNew = false;
  }
  for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
    e->isNew = false;
  }
  
  // Step A: Compute new positions from the original mesh.
  //  A.1 For each old vertex, compute its new position.
  for (VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++) {
    Vector3D neighborSum(0,0,0);
    int n = 0;
    // Iterate around v using its halfedge (neighbors given by h->twin()->vertex())
    HalfedgeIter h = v->halfedge();
    if (h == mesh.halfedgesEnd()) continue; // safety check
    HalfedgeIter start = h;
    do {
      VertexIter nb = h->twin()->vertex();
      neighborSum += nb->position;
      n++;
      h = h->twin()->next();
    } while (h != start);
    
    double u;
    if (n == 3)
      u = 3.0 / 16.0;
    else
      u = 3.0 / (8.0 * n);
    
    v->newPosition = (1 - n * u) * v->position + u * neighborSum;
  }
  
  //  A.2 For each original edge, compute the new midpoint position.
  // Process only original edges (isNew == false)
  for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
    if (e->isNew) continue;
    HalfedgeIter h = e->halfedge();
    VertexIter A = h->twin()->vertex();
    VertexIter B = h->vertex();
    VertexIter C = h->next()->vertex();
    VertexIter D = h->twin()->next()->vertex();
    e->newPosition = (3.0/8.0) * (A->position + B->position) +
                     (1.0/8.0) * (C->position + D->position);
  }
  
  // Step B: Subdivide the mesh.
  // Collect a list of original edges to split.
  vector<EdgeIter> originalEdges;
  for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
    if (!e->isNew)
      originalEdges.push_back(e);
  }
  
  // Split each original edge.
  // The splitEdge(e) function should insert a new vertex at the midpoint,
  // assign the new edge positions (from e->newPosition), and mark new elements as new.
  for (size_t i = 0; i < originalEdges.size(); i++) {
    EdgeIter e = originalEdges[i];
    VertexIter newV = mesh.splitEdge(e);
    if (newV != VertexIter()) {
      newV->isNew = true;
      // Also assign the new vertex position computed earlier from its parent edge.
      newV->newPosition = e->newPosition;
    }
  }
  
  // Step C: Flip new edges that connect an old vertex with a new vertex.
  // We iterate over all edges flagged as new.
  for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
    if (!e->isNew) continue;
    // Get endpoints.
    HalfedgeIter h = e->halfedge();
    VertexIter v1 = h->twin()->vertex(); // one endpoint
    VertexIter v2 = h->vertex();           // other endpoint
    // Flip the edge if one endpoint is old and the other is new.
    if ((v1->isNew && !v2->isNew) || (!v1->isNew && v2->isNew)) {
      mesh.flipEdge(e);
    }
  }
  
  // Step D: Update vertex positions.
  // Use the precomputed new positions for all vertices.
  for (VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++) {
    v->position = v->newPosition;
  }
}

} // namespace CGL  