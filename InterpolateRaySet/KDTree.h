#ifndef __KDTREE_H
#define __KDTREE_H

#include <array>
#include <vector>
#include <limits>
#include <stack>
#include <random>
#include <functional>
#include <iterator>
// #include <mutex>
// #include <atomic>

namespace KDTree
	{

	/* KDTree: 
	A KD tree is a tree structure, organizing a set of points in n-dimensional space. 
	For a given set of points (an unordered 1-dimensional array of n-dimensional points), creating the tree follows this procedure:
	1) The overall bounding box is determined
	2) The set of points is partitioned according to the coordinates of the first dimension, 
		such that we have two subsets of equal size (for an even number of points) or two subsets whose size differs by one;
		the "lower" half and the "upper" half.
	3) Each of the two subsets is partitioned into two halves, a lower half and an upper half, but now along the second dimension
	4) The procedure continues, cycling through the dimensions, until a stopping criterion is reached. In our case, we continue until 
		each subset contains exactly one point.
	The data structures are organized as follows:
	The tree structure consists of nodes. The full set of points is contained in the "root node", which has two child nodes: the upper and the lower half
	Each child node again has two child nodes, except if it contains only a single point. Nodes with single points and no child nodes are called "leaf nodes"
	
	We need a number of bookkeeping data structures.
	1) The underlying "point array" pts_, length m, of n-dim points, as a std::vector<std::array<float, ndim>>. Each array entry is a n-dim point.
	2) The same array transposed, as a std::array<std::vector<float>, ndim>, ptCoords_. Each array entry is a long vector (length m) 
		of a single coordinate. We use this transposed array for the partitioning algorithm.
	3) An array of integer indices into the "point array", idx_ of length m: the "point index array", which is a permutation of the range 0..(m-1). 
		The permutation is such that each node contains a consecutive range of entries in this "point index array". 
		This array essentially contains the sorting information in a way that the points belonging to any node have consecutive entries in this array.
		This allows to denote the potentially large number of points in a node to be stored by two numbers, 
		the first and the one-beyond-last point index entry.
		This array answers the question: "Which actual points belong to a node?"
	4) The reverse permutation of the "point index array", reverseIdx_: the "reverse point index array". Length = m.
		This array answers the question: "For a given point, with position i in the point array: Which entry in the point index array points to this i'th point?"
	5) The array of nodes, nodes_, starting with the root node, followed by the two child nodes of the root node, followed by the two 
		child nodes of the first child and so on.
		This "Node array" has 2 * m - 1 entries (draw the trees for m = 1, 2, 3, 4 to verify).
	6) An array of indices into the node array, nodeIdx_: The "node index array", of length m.
		This "node index array" answers the question: "Which node in the node array is the leaf node containing the i'th point in the point array?"
		It follows that i == point index array[node array[node index array[i]]->first point], or in other words:
		The i'th entry in the node index array is an integer k. The k'th entry in the node array is a node l. This node l has an integer data member, beginPoint. 
		The beginPoint'th entry in the point index array is an integer j. For a consistent tree, j == i.
	*/

	
	// Contains type definitions and the underlying dimension
	class Def
		{
		public:
			static constexpr size_t dim = 4; // max. dim == 8, restricted by the way we track edge nodes
			using TReal = float; // for ray files, float is appropriate -- all binary formats use float, all text formats don't typically use more than seven digits per number

			using TKDPoint = std::array<TReal, dim>;  // a single n-dim point
			using TKDPoints = std::vector<TKDPoint>;  // the type of the underlying array of points

			using TIdx = std::uint32_t; // ray files larger than 2^31 = 2 billion rays are not practical. If needed, change to size_t
			
			struct TPointIdx { TIdx pi_; }; // holds an index into the array of TKDPoint, pts_
			struct TIdxIdx { // holds an index into the array of TPointIdx, idx_: the next level of indirection. 
				TIdx ii_;	 // To be used as begin/end of range of elements within an array of TPointIdx.
				TIdxIdx& operator++() { ++ii_; return *this; };
				TIdxIdx operator++(int) { return TIdxIdx{ ii_++ }; };
				}; 
			struct TNodeIdx { TIdx ni_; }; // holds an index into the array of nodes, nodes_
			static constexpr TIdx invalidIdx = std::numeric_limits<TIdx>::max(); 
			using TPointIdxArray = std::vector<TPointIdx>; 
			using TIdxIdxArray = std::vector<TIdxIdx>;
			using TNodeIdxArray = std::vector<TNodeIdx>;
			using TRealArray = std::vector<TReal>;
			using TBox = std::pair<Def::TKDPoint, Def::TKDPoint>;
		};

	Def::TReal Distance(const Def::TKDPoint& p1, const Def::TKDPoint& p2);

	Def::TBox BoundingBox(const Def::TKDPoints& pts);

	bool IsInBox(const Def::TBox& boundingBox, const Def::TKDPoint& pt);
	bool AreInBox(const Def::TBox& boundingBox, const Def::TKDPoints& pts);
	std::array<Def::TKDPoint, 16> BoxCorners(const Def::TBox& box);
	bool BoxesOverlap(const Def::TBox& lhs, const Def::TBox& rhs); // touch is sufficient, i.e. neighboring boxes do overlap
	bool LhsBoxIsWithinRhsBox(const Def::TBox& lhs, const Def::TBox& rhs); // lhs is fully within closed rhs;  i.e. BoxIsWithinBox( b, b ) == true
	Def::TKDPoint MidPoint(const Def::TKDPoint& p0, const Def::TKDPoint& p1);


	struct TNode
		{
		using TReal = Def::TReal;
		// these two indices point into the idx_ member of TKDTree
		Def::TIdxIdx beginpt_; // [beginpt_;endpt_) is the range of indices of points which belong to this node
		Def::TIdxIdx endpt_; // if (endpt_ - beginpt_ == 1) then this node contains one point => leaf node
		// these two points are n-dim points with floating point coordinates
		Def::TKDPoint corner0_; // the two corners defining the cube containing all points of this node
		Def::TKDPoint corner1_;
		// these three indices point into the TNodeArray member of TKDTree
		// == invalidIdx == std::numeric_limits<TIdx>::max() in case there is no mother or child node
		Def::TNodeIdx motherNode_; // the mother node (this node is either low or hi child of the mother node)
		Def::TNodeIdx lowChildNode_; // the child node with lower values of current coordinate split
		Def::TNodeIdx hiChildNode_;// the child node with higher values of current coordinate split
		std::uint16_t splitDim_; // 0 .. Def::dim-1: The dimension this node is split into
		std::uint16_t isEdgeNode_; // bit field for edge faces: 
		// bit 0..(dim-1) is 1 if corner0_[0..dim] is boundary
		// bit 8 .. dim+7 is 1 if corner1_[dim] is boundary
		// all other bits. i.e. dim..7 and (dim+8)..15 are zero
		bool IsLeaf() const { return endpt_.ii_ - beginpt_.ii_ == 1; }
		Def::TIdx NPoints() const { return endpt_.ii_ - beginpt_.ii_; }
		Def::TReal Distance(const Def::TKDPoint& pt) const; // return distance of pt to nearest corner, 0 if inside
		Def::TReal Volume() const; // product over i of (corner1_[i] - corner0_[i]
		Def::TBox Box() const { return Def::TBox{ corner0_, corner1_ }; };

		// partition node into nPoints sub-blocks which all have same volume.
		// precondition: point must be within this node bounding box
		// returns list of points within sub-blocks. These points are all at center, except for the sub-block which contains pt
		// rv has nPoints-1 center points and pt
		std::vector<Def::TKDPoint> Partition(Def::TIdx nPoints, const Def::TKDPoint& pt) const;
		std::vector<Def::TKDPoint> RandomPartition(Def::TIdx nPoints, const Def::TKDPoint & pt, const std::function<TReal()>& ranGen) const;

		};

	class TKDTree: public Def
		{
		public:
			// constructor postcondition: pts_ is a copy of pts, idx_ is filled with 0,1,...,pts.size()-1, nodes_ is empty
			explicit TKDTree(const TKDPoints& pts);
			explicit TKDTree(TKDPoints&& pts);
			// create the tree. Precondition: freshly constructed TKDTree. Postcondition: Nodes array is filled, Index array is shuffled
			void CreateTree();
			void CheckConsistency() const; // check various invariants about internal state
			// reduce bounding boxes of leaf edge nodes such that edge is fac*(average node size) away from point
			void ShrinkEdgeNodes(TReal fac = 0.5);

			// read access functions
			const TKDPoints& Points() const; // the array of n-dim points
			using TNodeArray = std::vector<TNode>;
			const TNodeArray& Nodes() const; // the array of nodes which make up the tree
			const TPointIdxArray& Index() const; // the array of indices into the point array
			const TIdxIdxArray& ReverseIndex() const; // the inverse permutation of Index()
			const TNodeIdxArray& NodeIndex() const; // which leaf node contains the i'th point in Points()
			TBox BoundingBox() const;

			// make indices into various arrays type safe
			// convention: use only these access functions. Could be implemented by creating proxy classes for pts_,idx_,reverseIdx_ etc
			// but that is just too much effort.
			const TKDPoint& Point(TPointIdx pi) const { return pts_[pi.pi_]; };
			TKDPoint& Point(TPointIdx pi) { return pts_[pi.pi_]; };
			TPointIdx  Index(TIdxIdx ii) const { return idx_[ii.ii_]; };
			TPointIdx& Index(TIdxIdx ii) { return idx_[ii.ii_]; };
			TIdxIdx  ReverseIndex(TPointIdx pi) const { return reverseIdx_[pi.pi_]; };
			TIdxIdx& ReverseIndex(TPointIdx pi) { return reverseIdx_[pi.pi_]; };
			TNodeIdx  NodeIndex(TPointIdx pi) const { return nodeIdx_[pi.pi_]; }
			TNodeIdx& NodeIndex(TPointIdx pi) { return nodeIdx_[pi.pi_]; }


			const TNode& Node(TNodeIdx ni) const { return nodes_[ni.ni_]; };
			TNode& Node(TNodeIdx ni) { return nodes_[ni.ni_]; };


			// return the index of the node which contains pt. If out of bounds, returns invalidIdx. 
			// When exactly on border, returns node with lower coordinate values.
			TNodeIdx Locate(const TKDPoint& pt) const; 
			// returns the index of the node which contains point i; returns invalidIdx if i >= pts_.size()
			TNodeIdx Locate(TPointIdx pi) const;

			//// find the indices of the points which are within the box
			//TPointIdxArray LocatePointsWithinBox(const Def::TBox& box) const;
			// find the indices of the points which are within the box via KD tree
			TPointIdxArray LocatePointsWithinBox(const Def::TBox& box) const;

			// find the indices of the leaf nodes which overlap 
			TNodeIdxArray LocateOverlappingLeafNodes(const Def::TBox& box) const;


			// nearest neighbors
			struct TNearestNeighbors
				{
				explicit TNearestNeighbors(TIdx nNeighbors); // 
				TPointIdxArray i_points_; // the list of indices into pts_
				TNodeIdxArray i_nodes_; // the list of corresponding indices into nodes_
				TRealArray distances_; // the corresponding distances, sorted ascending
				};
			// for point pt, which may be anywhere (also outside), returns struct for n nearest neighbors
			TNearestNeighbors NearestNeighbors(const TKDPoint& pt, TIdx n) const;
			// for point i in the original points array
			TNearestNeighbors NearestNeighborsOfPoint(TPointIdx ipoint, TIdx n) const;
			// for node i in the node array, looks for neighbors of center point of node
			TNearestNeighbors NearestNeighborsOfNode(TNodeIdx inode, TIdx n) const;
			// for node i in the node array and point pt. precondition: pt is inside this node
			TNearestNeighbors NearestNeighborsOfNode(TNodeIdx inode, const TKDPoint& pt, TIdx n) const;

			TReal TotalVolume(const TNodeIdxArray& inodes) const;

		private:
			void Init();
			TNode RootNode();
			TNodeIdx AddNode(const TNode& n);
			void SetChildNodes(TNodeIdx mother, TNodeIdx loChild, TNodeIdx hiChild);
			void PushWork(TNodeIdx iNode);		
			TNodeIdx PopWork();
			void SplitNext();

			void LocatePointsOneNode(const TNode& node, const Def::TBox& box, Def::TPointIdxArray& pts) const;
			void LocateOverlappingNodesOneNode(const TNodeIdx& nodeIdx, const Def::TBox& box, Def::TNodeIdxArray& nodes) const;
			TKDPoints pts_;
			std::array<std::vector<Def::TReal>, Def::dim> ptCoords_; // rearranged pts_ for linear access of individual coordinates, avoids page faults in splitting algorithm
			Def::TBox boundingBox_;  // the overall bounding box
			TPointIdxArray idx_; // the shuffled index array which arranges pts into tree structure. Nodes contain contiguous ranges within this array
			TIdxIdxArray reverseIdx_; // the reverse shuffling: which idx_ entry points to the i'th point
			TNodeArray nodes_; // the list of nodes, starting with the root node. Contains 2 * pts_size() - 1 nodes
			TNodeIdxArray nodeIdx_; // which leaf node contains the i'th point in pts

			// used during construction, empty after construction
			std::stack<TNodeIdx> work_; // the nodes to be splitted
			// std::mutex mutex_;

			friend void TestKDTree4D();
		};

	void TestKDTree2D(std::string fn);

	void TestKDTree4D();



	} // namespace KDTree


#endif