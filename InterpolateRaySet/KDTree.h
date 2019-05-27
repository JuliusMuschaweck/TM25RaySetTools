#ifndef __KDTREE_H
#define __KDTREE_H

#include <array>
#include <vector>
#include <limits>
#include <stack>
// #include <mutex>
// #include <atomic>

namespace KDTree
	{

	class Def
		{
		public:
			using TReal = float;
			using TIdx = std::uint32_t;
			struct TPointIdx { TIdx pi_;}; // index into array of points
			struct TIdxIdx { TIdx ii_; }; // index into array of TPointIdx
			struct TNodeIdx { TIdx ni_; };
			static constexpr TIdx invalidIdx = std::numeric_limits<TIdx>::max();
			static constexpr size_t dim = 4;
			using TKDPoint = std::array<TReal, dim>;
			using TKDPoints = std::vector<TKDPoint>;
			using TIdxArray = std::vector<TIdx>;
			using TRealArray = std::vector<TReal>;
		};

	Def::TReal Distance(const Def::TKDPoint& p1, const Def::TKDPoint& p2);

	std::pair<Def::TKDPoint, Def::TKDPoint> BoundingBox(const Def::TKDPoints& pts);

	struct TNode
		{
		// these two indices point into the TNodeArray member of TKDTree
		Def::TIdx beginpt_; // [beginpt_;endpt_) is the range of indices of points which belong to this node
		Def::TIdx endpt_; // if (endpt_ - beginpt_ == 1) then this node contains one point => leaf node
		// these two points are n-dim points with floating point coordinates
		Def::TKDPoint corner0_; // the two corners defining the cube containing all points of this node
		Def::TKDPoint corner1_;
		// these three indices point into the TNodeArray member of TKDTree
		// == std::numeric_limits<TIdx>::max() in case there is no mother or child node
		Def::TIdx motherNode_; // the mother node (this node is either low or hi child of the mother node)
		Def::TIdx lowChildNode_; // the child node with lower values of current coordinate split
		Def::TIdx hiChildNode_;// the child node with higher values of current coordinate split
		std::uint16_t splitDim_; // 0 .. Def::dim-1: The dimension this node is split into
		std::uint16_t isEdgeNode_; // bit field for edge faces: 
		// bit 0..(dim-1) is 1 if corner0_[0..dim] is boundary
		// bit 8 .. dim+7 is 1 if corner1_[dim] is boundary
		// all other bits are zero
		bool IsLeaf() const { return endpt_ - beginpt_ == 1; }
		Def::TIdx NPoints() const { return endpt_ - beginpt_; }
		Def::TReal Distance(const Def::TKDPoint& pt) const; // return distance of pt to nearest corner, 0 if inside
		Def::TReal Volume() const;
		};

	class TKDTree: public Def
		{
		public:
			// constructor postcondition: pts_ is a copy of pts, idx_ is filled with 0,1,...,pts.size()-1, nodes_ is empty
			explicit TKDTree(const TKDPoints& pts);
			explicit TKDTree(TKDPoints&& pts);
			// create the tree. Precondition: freshly constructed TKDTree. Postcondition: Nodes array is filled, Index array is shuffled
			void CreateTree();
			// reduce bounding boxes of leaf edge nodes such that edge is fac*(average node size) away from point
			void ShrinkEdgeNodes(TReal fac = 0.5);

			// read access functions
			const TKDPoints& Points() const;
			using TNodeArray = std::vector<TNode>;
			const TNodeArray& Nodes() const;
			const TIdxArray& Index() const;
			const TIdxArray& ReverseIndex() const;
			std::pair<TKDPoint, TKDPoint> BoundingBox() const;


			// return the index of the node which contains pt. If out of bounds, returns invalidIdx. 
			// When exactly on border, returns node with lower coordinate values.
			TIdx Locate(const TKDPoint& pt) const; 
			// returns the index of the node which contains point i; returns invalidIdx if i >= pts_.size()
			TIdx Locate(TIdx i) const;

			// nearest neighbors
			struct TNearestNeighbors
				{
				explicit TNearestNeighbors(TIdx n); //
				TIdxArray i_points_; // the list of indices into pts_
				TIdxArray i_nodes_; // the list of corresponding indices into nodes_
				TRealArray distances_; // the corresponding distances, sorted ascending
				};
			// for point pt, which may be anywhere (also outside), returns struct for n nearest neighbors
			TNearestNeighbors NearestNeighbors(const TKDPoint& pt, TIdx n) const;
			// for point i in the original points array
			TNearestNeighbors NearestNeighborsOfPoint(TIdx ipoint, TIdx n) const;
			// for node i in the node array, looks for neighbors of center point of node
			TNearestNeighbors NearestNeighborsOfNode(TIdx inode, TIdx n) const;
			// for node i in the node array and point pt. precondition: pt is inside this node
			TNearestNeighbors NearestNeighborsOfNode(TIdx inode, const TKDPoint& pt, TIdx n) const;

			TReal TotalVolume(const TIdxArray& inodes) const;

		private:
			void Init();
			TNode RootNode();
			TIdx AddNode(const TNode& n);
			void SetChildNodes(TIdx mother, TIdx loChild, TIdx hiChild);
			void PushWork(TIdx iNode);		
			TIdx PopWork();
			void SplitNext();

			TKDPoints pts_;
			std::array<std::vector<Def::TReal>, Def::dim> ptCoords_; // rearranged pts_ for linear access of individual coordinates, avoids page faults in splitting algorithm
			std::pair<TKDPoint, TKDPoint> boundingBox_;
			TIdxArray idx_; // the shuffled index array which arranges pts into tree structure
			TIdxArray reverseIdx_; // the reverse shuffling: which idx_ entry points to the i'th point
			TNodeArray nodes_; // the list of nodes, starting with the root node. Contains 2 * pts_size() - 1 nodes
			TIdxArray nodeIdx_; // which leaf node contains the i'th point in pts
			std::stack<TIdx> work_; // the nodes to be splitted
			// std::mutex mutex_;
		};

	void TestKDTree2D(std::string fn);

	void TestKDTree4D();



	} // namespace KDTree


#endif