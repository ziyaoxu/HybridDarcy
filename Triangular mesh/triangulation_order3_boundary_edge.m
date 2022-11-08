function boundary_edge = triangulation_order3_boundary_edge ( triangle_node )

%*****************************************************************************80
%
%% TRIANGULATION_ORDER3_BOUNDARY_EDGE returns the boundary edges.
%
%  Discussion:
%
%    This routine is given a triangulation, an abstract list of triples
%    of nodes.  It is assumed that the nodes in each triangle are listed
%    in a counterclockwise order, although the routine should work 
%    if the nodes are consistently listed in a clockwise order as well.
%
%    It is assumed that each edge of the triangulation is either 
%    * an INTERIOR edge, which is listed twice, once with positive
%      orientation and once with negative orientation, or;
%    * a BOUNDARY edge, which will occur only once.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    12 January 2014
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer TRIANGLE_NUM, the number of triangles.
%
%    Input, integer TRIANGLE_NODE(3,TRIANGLE_NUM), the nodes that make up the
%    triangles.  These should be listed in counterclockwise order.
%
%    Output, integer BOUNDARY_EDGE(2,BOUNDARY_EDGE_NUM), the pairs of nodes
%    that make up the boundary edges.
%

%
%  Set up the edge array.
%
  triangle_num = size(triangle_node,2);
  edge_num = 3 * triangle_num;
  edge = zeros ( 3, edge_num );

  edge(1:2,               1:  triangle_num) = triangle_node(1:2,1:triangle_num);
  edge(1:2,  triangle_num+1:2*triangle_num) = triangle_node(2:3,1:triangle_num);
  edge(1  ,2*triangle_num+1:3*triangle_num) = triangle_node(3,  1:triangle_num);
  edge(2  ,2*triangle_num+1:3*triangle_num) = triangle_node(1,  1:triangle_num);
%
%  For sorting, we need to reorder some edges.  
%  But for recovery later, we need the original ordering.
%  So add a third item which is 1 or 0 depending on whether we reordered.
%
  edge(3,1:3*triangle_num) = ( edge(1,1:3*triangle_num) < edge(2,1:3*triangle_num) );
%
%  In each column, force the smaller entry to appear first.
%
  e1(1:edge_num) = min ( edge(1:2,1:edge_num) );
  e2(1:edge_num) = max ( edge(1:2,1:edge_num) );

  edge(1,1:edge_num) = e1(1:edge_num);
  edge(2,1:edge_num) = e2(1:edge_num);
%
%  Now ascending sort the column array.
%
  edge = ( sortrows ( edge' ) )';
%
%  The boundary edges are the elements that occur just once in EDGE.
%
  boundary_edge = zeros ( 2, edge_num );

  e = 0;
  be = 0;
  
  while ( 1 )
%
%  Only one edge left.
%
    if ( e == edge_num - 1 )
      e = e + 1;
      be = be + 1;
      if ( edge(3,e) )
        boundary_edge(1:2,be) = edge(1:2,e);
      else
        boundary_edge(1,be) = edge(2,e);
        boundary_edge(2,be) = edge(1,e);
      end
      break
    end
%
%  We can compare edges E+1 and E+2.
%
    if ( edge(1:2,e+1) == edge(1:2,e+2) )
      e = e + 2;
    else
      e = e + 1;
      be = be + 1;
      if ( edge(3,e) )
        boundary_edge(1:2,be) = edge(1:2,e);
      else
        boundary_edge(1,be) = edge(2,e);
        boundary_edge(2,be) = edge(1,e);
      end
    end
    
  end
%
%  Reduce the storage for BOUNDARY_EDGE.
%
  boundary_edge = boundary_edge(1:2,1:be);
  boundary_edge = int32(boundary_edge);
  return
end
