#include"KDTree.h"
#include<numeric>
#include<algorithm>
#include<limits>
#include<assert.h>
#include<random>
#include<fstream>
#include<chrono>
#include<Timer.h>
#include <iostream>


namespace KDTree
	{
	TKDTree::TKDTree(const TKDPoints& pts)
		: pts_(pts)
		{
		Init();
		}

	TKDTree::TKDTree(TKDPoints&& pts)
		: pts_(std::move(pts))
		{
		Init();
		}

	void TKDTree::CreateTree()
		{
		constexpr Def::TReal big = std::numeric_limits<Def::TReal>::max();
		// construct root node
//		std::cout << sizeof(TNode) << std::endl;
//		std::cout << "reserving " << 2 * pts_.size() * sizeof(TNode) << " bytes of memory for nodes " << std::endl;
		nodes_.reserve(2 * pts_.size()); // one more than needed
		TNode root = RootNode();
		TNodeIdx iRoot = AddNode(root);
		PushWork(iRoot);
		TTimer tim;
		tim.tic();
		int i = 0;
		double elapsed = 0;
		while (!work_.empty())
			{
			if ((++i) % 1000 == 0)
				{
				elapsed = tim.toc();
				if (elapsed > 1)
					{
					tim.tic();
					double frac = (1.0 * nodes_.size()) / ((2 * pts_.size() - 1));
					std::cout << nodes_.size() << "/" << (2 * pts_.size() - 1) << "(" << round(frac * 1000) / 10 << "%)";
					std::cout.flush();
					}
				}
			SplitNext();
			}
		TIdxIdx ii{ 0 };
		for (TPointIdx i : idx_)
			{
			ReverseIndex(i) = ii;
			++(ii.ii_);
			}
		std::cout << " done\n";
		}

	void TKDTree::CheckConsistency() const
		{
		auto Check = [](bool test, const char* message)
			{
			if (!test)
				throw std::runtime_error("TKDTree::CheckConsistency: " + std::string(message));
			};
		// uint32 can hold # of points
		Check(pts_.size() < std::numeric_limits<TIdx>::max(), "too many points");
		// ptCoords_ = transpose(pts_)
		size_t npts = pts_.size();
		for (size_t i = 0; i < Def::dim; ++i)
			{
			Check(npts == ptCoords_[i].size(), "ptCoords_ has wrong dimension");
			for (size_t j = 0; j < npts; ++j)
				Check(ptCoords_[i][j] == pts_[j][i], "ptCoords_ is not pts_ transposed");
			}
		// all pts within bounding box
		for (size_t j = 0; j < npts; ++j)
			Check(IsInBox(boundingBox_, pts_[j]), "point outside bounding box");
		Check(idx_.size() == npts, "idx_ has wrong size");
		// idx_ is a complete permutation of 0..n-1
		{
		std::vector<int>  ni(npts, 0);
		for (auto i : idx_)
			{
			Check(i.pi_ < npts, "idx_ contains illegal value");
			ni[i.pi_] += 1;
			}
		for (auto i : ni)
			Check(i == 1, "idx_ is not a complete permutation of 0..npts-1");
		}
		// reverseIdx_ is reverse permutation of idx_
		Check(reverseIdx_.size() == npts, "reverseIdx_ has wrong size");
		// idx_ is a complete permutation of 0..n-1
		{
		std::vector<int>  nri(npts, 0);
		for (auto ri : reverseIdx_)
			{
			Check(ri.ii_ < npts, "reverseIdx_ contains illegal value");
			Check(reverseIdx_[idx_[ri.ii_].pi_].ii_ == ri.ii_, "reverseIdx_ is not inverse permutation of idx_");
			}
		}
		// node structure forms proper tree:
		//		by traversing tree all nodes are visited exactly once
		//		childs have correct mother
		//		endpt_ - beginpt_ >= 1
		//		endpt_ < nPoints
		//		each point is contained in exactly one leaf node
		//		# of leaf nodes == # of points
		//		all points are within corner0_ and corner1_
		//		splitDim_ < dim
		//		isEdgeNode_ matches with child's 
		Check(nodes_.size() < std::numeric_limits<TIdx>::max(), "too many nodes");
		{
		std::vector<int>  nni(nodes_.size(), 0);
		std::vector<int>  nii(npts, 0);
		std::stack<TIdx> work;
		auto DoNode = [this, &nni, &nii, npts, &work, &Check](TIdx i)
			{
			nni[i] += 1;
			const TNode& node = this->Nodes()[i];
			if (node.lowChildNode_.ni_ != Def::invalidIdx)
				{
				work.push(node.lowChildNode_.ni_);
				Check(this->Nodes()[node.lowChildNode_.ni_].motherNode_.ni_ == i, "child node has wrong mother");
				}
			if (node.hiChildNode_.ni_ != Def::invalidIdx)
				{
				work.push(node.hiChildNode_.ni_);
				Check(this->Nodes()[node.hiChildNode_.ni_].motherNode_.ni_ == i, "child node has wrong mother");
				}
			Check(node.endpt_.ii_ > node.beginpt_.ii_, "node: endpt must be > beginpt");
			Check(node.endpt_.ii_ <= npts, "node: endpt > npts");
			if (node.IsLeaf())
				{
				nii[node.beginpt_.ii_] += 1;
				Check((node.lowChildNode_.ni_ == Def::invalidIdx) && (node.hiChildNode_.ni_ == Def::invalidIdx),
					"node: Leaf nodes must not have child nodes");
				}
			else
				{
				Check((node.lowChildNode_.ni_ != Def::invalidIdx) && (node.hiChildNode_.ni_ != Def::invalidIdx),
					"node: Non-Leaf nodes must have child nodes");
				}
			for (TIdx i = node.beginpt_.ii_; i < node.endpt_.ii_; ++i)
				{
				TKDPoint pt = pts_[Index()[i].pi_];
				for (int idim = 0; idim < Def::dim; ++idim)
					Check((node.corner0_[idim] <= pt[idim]) && (node.corner1_[idim] >= pt[idim]), "node: points must be in corner0_, corner1_ box");
				}
			Check(node.splitDim_ < Def::dim, "node: splitDim_ too large");

			};
		work.push(0);
		while (!work.empty())
			{
			TIdx next = work.top();
			work.pop();
			DoNode(next);
			}
		for (auto i : nni)
			Check(i == 1, "Nodes must appear exactly once in tree");
		for (auto i : nii)
			Check(i == 1, "Each point must appear exactly once in a leaf node");
		}
		// nodeIdx_ points to nodes which are leaf nodes and actually contain the point
		Check(NodeIndex().size() == npts, "nodeIndex_ must have npts entries");
		TPointIdx pi{ 0 };
		for (auto ni : NodeIndex())
			{
			const TNode& node = Node(ni);
			Check(node.IsLeaf(), "nodes in nodeIndex_ must be leaf nodes");
			Check(Index(node.beginpt_).pi_ == (pi.pi_)++, "node in noneIndex_ does not point to correct point");
			}
		// work_ is empty
		Check(work_.empty(), "work_ is not empty");
		}

	void TKDTree::ShrinkEdgeNodes(TReal fac)
		{
		Def::TReal avgPtsPerDim = pow(static_cast<Def::TReal>(pts_.size()), static_cast<Def::TReal>(1.0 / Def::dim));
		TKDPoint d;
		for (int i = 0; i < dim; ++i)
			d[i] = fac * (boundingBox_.second[i] - boundingBox_.first[i]) / avgPtsPerDim;
		for (auto& nd : nodes_)
			{
			if (nd.IsLeaf() && (nd.isEdgeNode_ != 0x0000))
				{
				TKDPoint pt = Point(Index(nd.beginpt_));
				std::uint16_t lo = 0x0001;
				std::uint16_t hi = 0x0100;
				for (int i = 0; i < dim; ++i)
					{
					if (nd.isEdgeNode_ & lo)
						nd.corner0_[i] = pt[i] - d[i];
					if (nd.isEdgeNode_ & hi)
						nd.corner1_[i] = pt[i] + d[i];
					lo = lo << 1;
					hi = hi << 1;
					}
				}
			}
		}

	const Def::TKDPoints& TKDTree::Points() const
		{
		return pts_;
		}

	const TKDTree::TNodeArray& TKDTree::Nodes() const
		{
		return nodes_;
		// TODO: insert return statement here
		}

	const Def::TPointIdxArray& TKDTree::Index() const
		{
		return idx_;
		}

	const Def::TIdxIdxArray& TKDTree::ReverseIndex() const
		{
		return reverseIdx_;
		}

	const TKDTree::TNodeIdxArray& TKDTree::NodeIndex() const
		{
		return nodeIdx_;
		}

	Def::TBox TKDTree::BoundingBox() const
		{
		return boundingBox_;
		}

	Def::TNodeIdx TKDTree::Locate(const TKDPoint& pt) const
		{
		const TNode* nd = &(nodes_[0]); // start with foot node
		if (!IsInBox(nd->Box(), pt))
			return TNodeIdx{ invalidIdx };
		TNodeIdx idx{ 0 };
		while (!(nd->IsLeaf()))
			{
			size_t sd = nd->splitDim_;
			if (pt[sd] < Node(nd->lowChildNode_).corner1_[sd])
				idx = nd->lowChildNode_;
			else
				idx = nd->hiChildNode_;
			nd = &(Node(idx));
			}
		return idx;
		}

	Def::TNodeIdx TKDTree::Locate(TPointIdx pi) const
		{ // could be simplified with introducing a "reverseNodeIdx_" data member
		if (pi.pi_ >= pts_.size())
			return TNodeIdx{ invalidIdx };
		TIdxIdx ri = ReverseIndex(pi);
		TNodeIdx nb{ 0 };
		while (!(Node(nb).IsLeaf()))
			{
			TNodeIdx dlo = Node(nb).lowChildNode_;
			if (ri.ii_ < Node(dlo).endpt_.ii_)
				nb = dlo;
			else
				nb = Node(nb).hiChildNode_;
			}
		return nb;
		}

	//Def::TPointIdxArray TKDTree::LocatePointsWithinBox(const Def::TBox& box) const
	//	{
	//	TPointIdxArray rv;
	//	size_t npts = Points().size();
	//	for (size_t i = 0; i < npts; ++i)
	//		{
	//		TPointIdx pi{ static_cast<Def::TIdx>(i) };
	//		if (IsInBox(box, Point(pi)))
	//			{
	//			rv.push_back(pi);
	//			}
	//		}
	//	return rv;
	//	}

	void TKDTree::LocatePointsOneNode(const TNode& node, const Def::TBox& box, Def::TPointIdxArray& pts) const
		{// calls itself recursively
		Def::TBox nodeBox{ node.corner0_, node.corner1_ };
		if (BoxesOverlap(box, nodeBox))
			{
 			if (LhsBoxIsWithinRhsBox(nodeBox, box)) // add all points within this box
				{
				for (TIdxIdx ii = node.beginpt_; ii.ii_ != node.endpt_.ii_; ++ii)
					{
					pts.push_back(Index(ii));
					}
				}
			else if (node.IsLeaf())
				{
				TIdxIdx ii = node.beginpt_;
				TPointIdx pi = Index(ii);
				const TKDPoint& p = Points()[pi.pi_];
				//if (IsInBox(box, p))
				//	{
				//	pts.push_back(pi);
				//	}
				pts.push_back(pi);
				}
			else
				{
				LocatePointsOneNode(Node(node.lowChildNode_), box, pts);
				LocatePointsOneNode(Node(node.hiChildNode_), box, pts);
				}
			}
		}

	void AddLeafNodes(const TKDTree& kdt, Def::TNodeIdx nodeIdx, Def::TNodeIdxArray& nodes)
		{
		const TNode& node = kdt.Node(nodeIdx);
		if (node.IsLeaf())
			nodes.push_back(nodeIdx);
		else
			{
			AddLeafNodes(kdt, node.lowChildNode_, nodes);
			AddLeafNodes(kdt, node.hiChildNode_, nodes);
			}
		}

	void TKDTree::LocateOverlappingNodesOneNode(const TNodeIdx& nodeIdx, const Def::TBox& box, Def::TNodeIdxArray& nodes) const
		{// calls itself recursively
		const TNode& node = Node(nodeIdx);
		Def::TBox nodeBox = node.Box();
		if (BoxesOverlap(box, nodeBox))
			{
			if (node.IsLeaf())
				{
				nodes.push_back(nodeIdx);
				}
			else if (LhsBoxIsWithinRhsBox(nodeBox, box)) // add all points within this box
				{
				AddLeafNodes(*this, nodeIdx, nodes);
				}
			else
				{
				LocateOverlappingNodesOneNode(node.lowChildNode_, box, nodes);
				LocateOverlappingNodesOneNode(node.hiChildNode_, box, nodes);
				}
			}
		}


	Def::TPointIdxArray TKDTree::LocatePointsWithinBox(const Def::TBox& box) const
		{
		TPointIdxArray rv;
		const TNode& root = Nodes().front();
		LocatePointsOneNode(root, box, rv);
		return rv;
		}


	Def::TNodeIdxArray TKDTree::LocateOverlappingLeafNodes(const Def::TBox& box) const
		{
		TNodeIdxArray rv;
		LocateOverlappingNodesOneNode(TNodeIdx{ 0 }, box, rv);
		return rv;
		}

	// return the index of the nearest node with >= n points if pt is outside bounding box of all nodes
	Def::TNodeIdx ClosestNodeIndex(const TKDTree::TNodeArray& nodes, const Def::TKDPoint& pt, Def::TIdx n)
		{
		Def::TNodeIdx inode{ 0 };
		for (;;)
			{
			const TNode& thisnode = nodes[inode.ni_];
			const TNode& lonode = nodes[thisnode.lowChildNode_.ni_];
			const TNode& hinode = nodes[thisnode.hiChildNode_.ni_];
			if (lonode.NPoints() < n && hinode.NPoints() < n)
				break;
			Def::TReal distlo = lonode.Distance(pt);
			Def::TReal disthi = hinode.Distance(pt);
			if (distlo < disthi)
				{
				if (lonode.NPoints() < n)
					break;
				else
					{
					inode = thisnode.lowChildNode_;
					if (lonode.NPoints() == n)
						break;
					continue;
					}
				}
			else
				{
				if (hinode.NPoints() < n)
					break;
				else
					{
					inode = thisnode.hiChildNode_;
					if (hinode.NPoints() == n)
						break;
					continue;
					}
				}
			}
		return inode;
		}

	void Heapify(TKDTree::TNearestNeighbors& nn, Def::TIdx n)
		{ // assume only element at 0 is wrong. Consider only first n elements
		Def::TIdx j0 = 0;
		Def::TIdx jTest = 1; //
		while (jTest < n)
			{
			if (jTest < (n - 1) && nn.distances_[jTest] < nn.distances_[jTest + 1])
				++jTest; // right underling is larger
			if (nn.distances_[j0] >= nn.distances_[jTest])
				break; // top is not smaller than largest underling -> stop
			// demote top -> swap top with underling
			std::swap(nn.distances_[j0], nn.distances_[jTest]);
			std::swap(nn.i_nodes_[j0], nn.i_nodes_[jTest]);
			std::swap(nn.i_points_[j0], nn.i_points_[jTest]);
			j0 = jTest;
			jTest = 2 * jTest + 1; // next underling
			}
		}

	TKDTree::TNearestNeighbors TKDTree::NearestNeighbors(const TKDPoint& pt, TIdx n) const
		{
		if (n >= pts_.size())
			throw std::runtime_error("TKDTree::NearestNeighbors: too many points requested");
		TNodeIdx inode = Locate(pt);
		if (inode.ni_ == invalidIdx)
			{ // pt is outside the bounding box
			// look for closest box with >= n points
			inode = ClosestNodeIndex(nodes_, pt, n);
			}
		return NearestNeighborsOfNode(inode, pt, n);
		}

	TKDTree::TNearestNeighbors TKDTree::NearestNeighborsOfPoint(TPointIdx ipoint, TIdx n) const
		{
		if (ipoint.pi_ >= pts_.size())
			throw std::runtime_error("TKDTree::NearestNeighborsOfPoint: ipoint > npoints");
		TNodeIdx inode = NodeIndex(ipoint);
		return NearestNeighborsOfNode(inode, pts_[ipoint.pi_], n);
		}

	TKDTree::TNearestNeighbors TKDTree::NearestNeighborsOfNode(TNodeIdx inode, TIdx n) const
		{
		// use center point of node
		const TNode& node = Node(inode);
		const TKDPoint& lo = node.corner0_;
		const TKDPoint& hi = node.corner1_;
		TKDPoint pt;
		for (int i = 0; i < Def::dim; ++i)
			pt[i] = static_cast<TReal>(0.5) * (lo[i] + hi[i]);
		return NearestNeighborsOfNode(inode, pt, n);
		}

	TKDTree::TNearestNeighbors TKDTree::NearestNeighborsOfNode(TNodeIdx inode, const TKDPoint& pt, TIdx n) const
		{
		if (inode.ni_ == invalidIdx)
			throw std::runtime_error("TKDTree::NearestNeighborsOfNode: invalidIdx");
		if (inode.ni_ >= nodes_.size())
			throw std::runtime_error("TKDTree::NearestNeighborsOfNode: inode out of range");
		// find smallest mother box with n points
		while (Node(inode).NPoints() < n)
			{
			if (inode.ni_ != 0)
				inode = Node(inode).motherNode_;
			else // root node has less than n points
				{
				n = Node(inode).NPoints();
				break;
				}
			}
		// now inode points to nearest box with >= n points
		// select n nearest within this, as heap
		const TNode& nd = Node(inode);
		TNearestNeighbors rv(n); // filled with BIG and invalid indices
		for (TIdxIdx ipt = nd.beginpt_; ipt.ii_ < nd.endpt_.ii_; ++(ipt.ii_))
			{
			const TKDPoint& thispt = Point(Index(ipt));
			TReal d = Distance(thispt, pt);
			if (d < rv.distances_[0])
				{
				// replace top of heap with new values
				rv.distances_[0] = d;
				rv.i_points_[0] = Index(ipt);
				rv.i_nodes_[0] = NodeIndex(Index(ipt));
				if (n > 1)
					Heapify(rv, n); // maintain heap structure
				}
			}
		// rv is filled with heap of n candidate points, all from same box
		// now traverse tree, considering only potentially better boxes
		std::stack<TNodeIdx> task;
		task.push(TNodeIdx{ 0 }); // start working with top node
		while (!task.empty())
			{
			TNodeIdx itodo = task.top();
			task.pop();
			if (itodo.ni_ == inode.ni_) // this whole node is taken care of already
				continue;
			const TNode& nTodo = Node(itodo);
			if (nTodo.Distance(pt) < rv.distances_[0]) // skip whole node if too far away
				{ // candidate node
				if (nTodo.IsLeaf())
					{ // consider point
					TPointIdx ipt = Index(nTodo.beginpt_);
					TReal d = Distance(Point(ipt), pt);
					if (d < rv.distances_[0])
						{// replace top of heap with new values
						rv.distances_[0] = d;
						rv.i_points_[0] = ipt;
						rv.i_nodes_[0] = itodo;
						if (n > 1)
							Heapify(rv, n); // maintain heap structure
						}
					}
				else
					{// put childs on task list
					task.push(nTodo.lowChildNode_);
					task.push(nTodo.hiChildNode_);
					}
				}
			}
		// now the heap of return values is filled
		// sort it
		TIdx nRemain = n;
		while (nRemain > 1)
			{
			--nRemain;
			std::swap(rv.distances_[0], rv.distances_[nRemain]);
			std::swap(rv.i_nodes_[0], rv.i_nodes_[nRemain]);
			std::swap(rv.i_points_[0], rv.i_points_[nRemain]);
			Heapify(rv, nRemain);
			}
		return rv;
		}

	Def::TReal TKDTree::TotalVolume(const TNodeIdxArray& inodes) const
		{
		TReal rv = 0;
		for (TNodeIdx i : inodes)
			{
			rv += Node(i).Volume();
			}
		return rv;
		}

	void TKDTree::Init()
		{
		idx_.resize(pts_.size());
		// std::iota(idx_.begin(), idx_.end(), 0);
		TIdx iota = 0;
		for (TPointIdx& pi : idx_)
			{
			pi.pi_ = iota++;
			}
		//reverseIdx_ = idx_;
		reverseIdx_.resize(pts_.size());
		// std::iota(idx_.begin(), idx_.end(), 0);
		iota = 0;
		for (TIdxIdx& ii : reverseIdx_)
			{
			ii.ii_ = iota++;
			}
		nodeIdx_.resize(pts_.size());
		std::fill(nodeIdx_.begin(), nodeIdx_.end(), TNodeIdx{ invalidIdx });
		for (size_t d = 0; d < Def::dim; ++d)
			ptCoords_[d].reserve(pts_.size());
		constexpr Def::TReal big = std::numeric_limits<Def::TReal>::max();
		Def::TKDPoint p0; p0.fill(big);
		Def::TKDPoint p1; p1.fill(-big);
		for (const auto& p : pts_)
			{
			for (size_t d = 0; d < Def::dim; ++d)
				{
				TReal pd = p[d];
				ptCoords_[d].push_back(pd);
				if (p0[d] > pd) p0[d] = pd;
				if (p1[d] < pd) p1[d] = pd;
				}
			}
		boundingBox_ = std::make_pair(p0, p1);
		}

	TNode TKDTree::RootNode()
		{
		TNode root{};
		root.beginpt_ = TIdxIdx{ 0 };
		root.endpt_ = TIdxIdx{ static_cast<TIdx>(pts_.size()) };
		std::pair<TKDTree::TKDPoint, TKDTree::TKDPoint> bb = KDTree::BoundingBox(pts_);
		root.corner0_ = bb.first; 
		root.corner1_ = bb.second;
		Def::TReal avgPtsPerDim = pow(static_cast<Def::TReal>(pts_.size()), static_cast<Def::TReal>(1.0 / Def::dim));
		for (int i = 0; i < Def::dim; ++i)
			{
			Def::TReal d = (root.corner1_[i] - root.corner0_[i]) / avgPtsPerDim * static_cast<Def::TReal>(0.1);
			root.corner0_[i] -= d;
			root.corner1_[i] += d;
			}
		root.motherNode_ = TNodeIdx{ Def::invalidIdx }; // root node is the only node w/o mother
		root.lowChildNode_ = TNodeIdx{ Def::invalidIdx }; // will be created later
		root.hiChildNode_ = TNodeIdx{ Def::invalidIdx }; // will be created later
		root.splitDim_ = 0; // we start with splitting along dimension 1
		std::array<uint16_t, 9> edgeFlags = { 0x0000, 0x0101, 0x0303, 0x0707, 0x0f0f, 0x1f1f, 0x3f3f, 0x7f7f, 0xffff };
		root.isEdgeNode_ = edgeFlags[dim];
		return root;
		}

	Def::TNodeIdx TKDTree::AddNode(const TNode& n)
		{
		// std::lock_guard<std::mutex> lock(mutex_);
		// this is where thread safe node adding may come in, e.g. by locking a mutex
		nodes_.push_back(n);
		return TNodeIdx{ static_cast<TIdx>(nodes_.size()) - 1 };
		}

	void TKDTree::SetChildNodes(TNodeIdx mother, TNodeIdx loChild, TNodeIdx hiChild)
		{
		Node(mother).lowChildNode_ = loChild;
		Node(mother).hiChildNode_ = hiChild;
		}

	void TKDTree::PushWork(TNodeIdx iNode)
		{
		work_.push(iNode);
		}

	Def::TNodeIdx TKDTree::PopWork()
		{
		Def::TNodeIdx rv{ work_.top() };
		work_.pop();
		return rv;
		}

	void TKDTree::SplitNext()
		{
		// retrieve node to work on
		TNodeIdx next = PopWork();
		TNode node = Node(next);
		// the coordinates according to which the split shall occur
		// #ifndef NDEBUG
		// std::cout << "working on node " << next << '\n';
		// #endif
		// the 1-dim array of node.splitDim_'th point coordinates
		// we sort according to the values in this array.
		const std::vector<Def::TReal>& coords = ptCoords_[node.splitDim_]; 
		// the range of indices to be split
		TIdx nPts = node.endpt_.ii_ - node.beginpt_.ii_;
		assert(nPts >= 2);
		auto ibegin = idx_.begin() + node.beginpt_.ii_;
		auto iend = idx_.begin() + node.endpt_.ii_;
		auto imid = ibegin + nPts / 2;
		// the comparison function for nth_element
		auto comp = [&coords](TPointIdx i1, TPointIdx i2) -> bool
			{
			// return coords[this->idx_[i1]] < coords[this->idx_[i2]];
			return coords[i1.pi_] < coords[i2.pi_];
			};
		// the actual splitting: the range [ibegin, iend) of idx_ is rearranged
		std::nth_element(ibegin, imid, iend, comp);
		// now  [ibegin, imid) has low values of coords, and [imid, iend) high values.
		// imid points to the split point. But we want to split the bounding box between
		// the split point and its nearest lower neighbor
		auto closestLoIdx = std::max_element(ibegin, imid, comp);
		TReal closestLoCoord = coords[(*closestLoIdx).pi_];
		TReal splitPointCoord = coords[(*imid).pi_];
		TReal BBSplitCoord = (closestLoCoord + splitPointCoord) / 2;
		TNode dlo;
		TNode dhi;
		dlo.beginpt_ = node.beginpt_;
		dlo.endpt_.ii_ = node.beginpt_.ii_ + nPts / 2;
		dhi.beginpt_ = dlo.endpt_;
		dhi.endpt_ = node.endpt_;
		dlo.corner0_ = dhi.corner0_ = node.corner0_;
		dlo.corner1_ = dhi.corner1_ = node.corner1_;
		dlo.corner1_[node.splitDim_] = dhi.corner0_[node.splitDim_] = BBSplitCoord;
		dlo.motherNode_ = dhi.motherNode_ = next;
		dlo.lowChildNode_ = dhi.lowChildNode_ = dlo.hiChildNode_ = dhi.hiChildNode_ = TNodeIdx{ Def::invalidIdx };
		dlo.splitDim_ = dhi.splitDim_ = (node.splitDim_ + 1) % dim;
		dlo.isEdgeNode_ = dhi.isEdgeNode_ = node.isEdgeNode_; // initialize the edge node field of childs to this one
		// adjust bit field of dim
		// dlo.corner1_[dim] is not an edge any more
		std::uint16_t loFlag = 1;
		loFlag = loFlag << (8 + node.splitDim_);
		loFlag = ~loFlag;
		dlo.isEdgeNode_ &= loFlag;
		// dhi.corner0_[dim] is not an edge any more
		std::uint16_t hiFlag = 1;
		hiFlag = hiFlag << node.splitDim_;
		hiFlag = ~hiFlag;
		dhi.isEdgeNode_ &= hiFlag;
		TNodeIdx ilo = AddNode(dlo);
		TNodeIdx ihi = AddNode(dhi);
		SetChildNodes(next, ilo, ihi);
		if (!dlo.IsLeaf())
			PushWork(ilo);
		else
			NodeIndex(Index(dlo.beginpt_)) = ilo;
		if (!dhi.IsLeaf())
			PushWork(ihi);
		else
			NodeIndex(Index(dhi.beginpt_)) = ihi;
		}

	Def::TReal Distance(const Def::TKDPoint& p1, const Def::TKDPoint& p2)
		{
		auto sqr = [](Def::TReal r) {return r * r; };
		Def::TReal rv = 0;
		for (int i = 0; i < p1.size(); ++i)
			{
			rv += sqr(p2[i] - p1[i]);
			}
		return sqrt(rv);
		}

	std::pair<Def::TKDPoint, Def::TKDPoint> BoundingBox(const Def::TKDPoints& pts)
		{
		constexpr Def::TReal big = std::numeric_limits<Def::TReal>::max();
		Def::TKDPoint p0; p0.fill(big);
		Def::TKDPoint p1; p1.fill(-big);
		for (const auto& p : pts)
			{
			for (size_t d = 0; d < Def::dim; ++d)
				{
				if (p0[d] > p[d]) p0[d] = p[d];
				if (p1[d] < p[d]) p1[d] = p[d];
				}
			}
		return std::make_pair(p0, p1);
		}

	bool IsInBox(const Def::TBox& boundingBox, const Def::TKDPoint& pt)
		{
		for (int idim = 0; idim < pt.size(); ++idim)
			{
			if (pt[idim] < boundingBox.first[idim]) return false;
			if (pt[idim] > boundingBox.second[idim]) return false;
			}
		return true;
		}

	bool AreInBox(const Def::TBox& boundingBox, const Def::TKDPoints& pts)
		{
		for (auto& p : pts)
			{
			if (!IsInBox(boundingBox, p)) return false;
			}
		return true;
		}

	std::array<Def::TKDPoint, 16> BoxCorners(const Def::TBox& box)
		{
		const Def::TKDPoint& b0 = box.first;
		const Def::TKDPoint& b1 = box.second;
		std::array<Def::TKDPoint, 16> rv
			{
			Def::TKDPoint{b0[0], b0[1], b0[2], b0[3]},
			Def::TKDPoint{b0[0], b0[1], b0[2], b1[3]},
			Def::TKDPoint{b0[0], b0[1], b1[2], b0[3]},
			Def::TKDPoint{b0[0], b0[1], b1[2], b1[3]},
			Def::TKDPoint{b0[0], b1[1], b0[2], b0[3]},
			Def::TKDPoint{b0[0], b1[1], b0[2], b1[3]},
			Def::TKDPoint{b0[0], b1[1], b1[2], b0[3]},
			Def::TKDPoint{b0[0], b1[1], b1[2], b1[3]},
			Def::TKDPoint{b1[0], b0[1], b0[2], b0[3]},
			Def::TKDPoint{b1[0], b0[1], b0[2], b1[3]},
			Def::TKDPoint{b1[0], b0[1], b1[2], b0[3]},
			Def::TKDPoint{b1[0], b0[1], b1[2], b1[3]},
			Def::TKDPoint{b1[0], b1[1], b0[2], b0[3]},
			Def::TKDPoint{b1[0], b1[1], b0[2], b1[3]},
			Def::TKDPoint{b1[0], b1[1], b1[2], b0[3]},
			Def::TKDPoint{b1[0], b1[1], b1[2], b1[3]}
			};
		return rv;
		}

	bool BoxesOverlap(const Def::TBox& lhs, const Def::TBox& rhs)
		{
		for (size_t i = 0; i < 4; ++i)
			{
			if (lhs.first[i] >= rhs.second[i]
				|| lhs.second[i] <= rhs.first[i])
				{
				return false;
				}
			}
		return true;
		}

	bool LhsBoxIsWithinRhsBox(const Def::TBox& lhs, const Def::TBox& rhs)
		{// lhs is fully within closed rhs; 
		for (size_t i = 0; i < Def::dim; ++i)
			{
			if (lhs.first[i] <= rhs.first[i]
				|| lhs.second[i] >= rhs.second[i])
				{
				return false;
				}
			}
		return true;
		}

	Def::TKDPoint MidPoint(const Def::TKDPoint& p0, const Def::TKDPoint& p1)
		{
		Def::TReal _5 = static_cast<Def::TReal>(0.5);
		return Def::TKDPoint{ _5 * (p0[0] + p1[0]),_5 * (p0[1] + p1[1]), _5 * (p0[2] + p1[2]), _5 * (p0[3] + p1[3]) };
		}


	void TestKDTree2D(std::string fn)
		{
		assert(Def::dim == 2);
		size_t npts = 1000;
		Def::TKDPoints pts;
		std::mt19937 gen;
		gen.seed();
		std::uniform_real_distribution<float> dis(-1.0f, 1.0f);
		for (size_t i = 0; i < npts; ++i)
			{
			auto pt = Def::TKDPoint{ 2 * dis(gen), dis(gen) };
			if ((pt[0] * pt[0]) / 4 + pt[1] * pt[1] < 1)
				pts.push_back(pt);
			}

		std::mt19937 gen2;
		gen2.seed();

		std::vector<float> arr;
		arr.reserve(npts);
		for (int i = 0; i < npts; ++i)
			arr.push_back(dis(gen2));
		std::chrono::high_resolution_clock clock;
		auto tic = clock.now();
		auto toc = clock.now();
		auto dt = toc - tic;

		npts = pts.size();

		tic = clock.now();
		std::nth_element(arr.begin(), arr.begin() + npts / 2, arr.end());
		toc = clock.now();
		dt = toc - tic;
		std::cout << "elapsed " << std::chrono::duration_cast<std::chrono::microseconds>(dt).count() * 1e-6 << " seconds for nth_element with " << npts << " points\n";

		TKDTree kdt(pts);
		tic = clock.now();
		kdt.CreateTree();
		toc = clock.now();
		dt = toc - tic;
		std::cout << "elapsed " << std::chrono::duration_cast<std::chrono::microseconds>(dt).count() * 1e-6 << " seconds for KD tree generation with " << npts << " points\n";

		tic = clock.now();
		kdt.ShrinkEdgeNodes();
		toc = clock.now();
		dt = toc - tic;
		std::cout << "elapsed " << std::chrono::duration_cast<std::chrono::microseconds>(dt).count() * 1e-6 << " seconds for shrinking edge nodes of " << npts << " points\n";

		Def::TRealArray vols(npts);
		tic = clock.now();
		for (int i = 0; i < npts; ++i)
			{
			auto nn = kdt.NearestNeighbors(pts[i], 10);
			vols[i] = kdt.TotalVolume(nn.i_nodes_);
			}
		toc = clock.now();
		dt = toc - tic;
		std::cout << "elapsed " << std::chrono::duration_cast<std::chrono::microseconds>(dt).count() * 1e-6 << " seconds for volume of 10 neighbor nodes of " << npts << " points\n";

		const TKDTree::TNodeArray& nd = kdt.Nodes();

		// return; 
		std::ofstream f(fn);
		f << "% output of TestKDTree2D\n";
		f << "clear; figure(1); clf; hold on;\n";
		for ( auto const& n : nd)
			{
			if (n.IsLeaf())
				{
				Def::TKDPoint pt = pts[(kdt.Index()[n.beginpt_.ii_]).pi_];
				Def::TKDPoint bblo = n.corner0_;
				Def::TKDPoint bbhi = n.corner1_;
				f << "plot([" << bblo[0] << ',' << bbhi[0] << ',' << bbhi[0] << ',' << bblo[0] << ',' << bblo[0] << "], ";
				f << '[' << bblo[1] << ',' << bblo[1] << ',' << bbhi[1] << ',' << bbhi[1] << ',' << bblo[1] << "],'k');\n";
				f << "scatter(" << pt[0] << "," << pt[1] << ",'kx');\n";
				}
			}
		for (auto const& n : nd)
			{
			if (n.IsLeaf())
				{
				Def::TKDPoint pt = pts[(kdt.Index()[n.beginpt_.ii_]).pi_];
				Def::TKDPoint bblo = n.corner0_;
				Def::TKDPoint bbhi = n.corner1_;
				if (n.isEdgeNode_ & 0x0001) // lo x
					f << "plot([" << bblo[0] << ',' << bblo[0] << "],[" << bblo[1] << ',' << bbhi[1] << "],'r');\n";
				if (n.isEdgeNode_ & 0x0002) // lo y
					f << "plot([" << bblo[0] << ',' << bbhi[0] << "],[" << bblo[1] << ',' << bblo[1] << "],'g');\n";
				if (n.isEdgeNode_ & 0x0100) // hi x
					f << "plot([" << bbhi[0] << ',' << bbhi[0] << "],[" << bblo[1] << ',' << bbhi[1] << "],'b');\n";
				if (n.isEdgeNode_ & 0x0200) // hi y
					f << "plot([" << bblo[0] << ',' << bbhi[0] << "],[" << bbhi[1] << ',' << bbhi[1] << "],'m');\n";
				}
			}
		int nn = 10;
		TKDTree::TNearestNeighbors nn00 = kdt.NearestNeighbors({ 0,0 }, nn);
		TKDTree::TNearestNeighbors nn22 = kdt.NearestNeighbors({ 2.1f,1.1f }, nn);
		f << "scatter(0,0,'bs');\n";
		f << "scatter(2.1,1.1,'bs');\n";
		for (int i = 0; i < nn; ++i)
			{
			auto pt1 = pts[nn00.i_points_[i].pi_];
			auto pt2 = pts[nn22.i_points_[i].pi_];
			f << "scatter([" << pt1[0] << ',' << pt2[0] << "],[" << pt1[1] << ',' << pt2[1] << "],'m');\n";
			}
		f << "axis equal;\n";
		}


	void TestKDTree4D()
		{
		/* postcondition of CreateTree :
		* 1) pts_, ptCoords_[0..3], idx_, reverseIdx_, nodeIdx_ all have same length m
		* 2) nodes_ has size 2 * m - 1
		* 3) work_ is empty
		* 4) min/max of all pts_ are within boundingBox_ limits, which are somewhat larger
		* 4) idx_ is a full permutation of 0..m-1
		* 5) reverseIdx_[idx_[i]] == i for all i in 0..m-1
		* 6) nodes_ contains exactly m leaf nodes (with one point)
		* 7) all leaf nodes have no child nodes
		* 9) root node is all edges
		* for all non-leaf nodes:
		*	10) test correct splitting of bounding box into child nodes
		*	lowChildNode.corner1_[splitDim_] == hiChildNode.corner0_[splitDim_]
		*	lowChildNode.corner0_[splitDim_] == corner0_[splitDim_]
		*	hiChildNode.corner1_[splitDim_] == corner1_[splitDim_]
		*	for j = 0..dim-1 except splitDim:
		*		lowChildNode.corner0_[j] == hiChildNode.corner0_[j] == corner0_[j]
		*		lowChildNode.corner1_[j] == hiChildNode.corner1_[j] == corner1_[j]
		*	11) test correct edge node treatment
		*	for j = 0..dim-1 except splitDim:
		*		isEdgeNode_ values same for this node and both child nodes
		*	for j = splitDim
		*		isEdgeNode_ values false for low child corner1_ and high child corner0_
		* for all nodes:
		*	12) all points are within bounding box
		* for i = 0..m-1
		*	13) test if nodeIdx_ entry i really points to the leaf node which contains the single point index i
		*		set inode = Node(NodeIndex(TPointIdx{i}))
		*		then inode must be leaf node
		*		and inode.beginpt_.pi_ == i
		*/
		auto check = [](bool condition, const std::string& msg)
			{
			if (!condition)
				throw std::logic_error("TestKDTree4D: condition failed: " + msg);
			};

		// prepare
		check(Def::dim == 4,"dim == 4");
		size_t npts = 100;
		Def::TKDPoints pts;
		std::mt19937 gen;
		gen.seed();
		std::uniform_real_distribution<float> dis(-1.0f, 1.0f);
		for (size_t i = 0; i < npts; ++i)
			{
			Def::TKDPoint pt;
			std::generate(pt.begin(), pt.end(), [&]() {return dis(gen); });
			if (pt[0] * pt[0] + pt[1] * pt[1] < 1) // within unit circle
				pts.push_back(pt);
			}

		std::chrono::high_resolution_clock clock;
		auto tic = clock.now();
		auto toc = clock.now();
		auto dt = toc - tic;

		TKDTree kdt(pts);
		tic = clock.now();
		kdt.CreateTree();
		toc = clock.now();
		dt = toc - tic;


		// test
		// 1) pts_, ptCoords_[0..3], idx_, reverseIdx_, nodeIdx_ all have same length m
		size_t m = kdt.pts_.size();
		for (size_t i = 0; i < 4; ++i) check(m == kdt.ptCoords_[i].size(), "ptCoords_[i].size() == m");
		check(m == kdt.idx_.size(), "idx_.size() == m");
		check(m == kdt.reverseIdx_.size(), "reverseIdx_.size() == m");
		check(m == kdt.nodeIdx_.size(), "nodeIdx_.size() == m");
		// 2) nodes_ has size 2 * m - 1
		check(2*m-1 == kdt.nodes_.size(), "nodes_.size() == 2*m-1");
		// 3) work_ is empty
		check(kdt.work_.empty(), "work_ is empty");
		// 4) min / max of all pts_ are within respective boundingBox_ limits, which are somewhat larger
		auto actualBB = BoundingBox(kdt.pts_);
		for (size_t i = 0; i < Def::dim; ++i)
			{
			check(kdt.boundingBox_.first[i] <= actualBB.first[i],"bounding box lower limit of dim" + std::to_string(i));
			check(kdt.boundingBox_.second[i] >= actualBB.second[i], "bounding box upper limit of dim" + std::to_string(i));
			}
		// 4) idx_ is a full permutation of 0..m - 1
		std::vector<bool> test;
		test.resize(m);
		std::fill(test.begin(), test.end(), false);
		for (size_t i = 0; i < m; ++i)
			test[kdt.idx_[i].pi_] = true;
		check(std::all_of(test.begin(), test.end(), [](bool b) {return b; }),"idx_ is a full permutation of 0..m - 1");
		// 5) reverseIdx_[idx_[i]] == i for all i in 0..m - 1
		std::fill(test.begin(), test.end(), false);
		for (size_t i = 0; i < m; ++i)
			test[i] = (kdt.reverseIdx_[kdt.idx_[i].pi_].ii_ == i);
		check(std::all_of(test.begin(), test.end(), [](bool b) {return b; }), "reverseIdx_[idx_[i]] == i for all i in 0..m - 1");
		// 6) nodes_ contains exactly m leaf nodes(with one point)
		size_t nLeaf = std::count_if(kdt.nodes_.begin(), kdt.nodes_.end(), [](const TNode& node) {return node.IsLeaf(); });
		check(nLeaf == m, "nodes_ contains exactly m leaf nodes(with one point)");

		bool rootNodeFound = false;
		for (const auto& node : kdt.nodes_)
			{
			// 7) all leaf nodes have no child nodes
			if (node.IsLeaf())
				{
				check(node.lowChildNode_.ni_ == Def::invalidIdx, "all leaf nodes have no low child node");
				check(node.hiChildNode_.ni_ == Def::invalidIdx, "all leaf nodes have no hi child node");
				}
			//8) only nodes_[0] has no mother node
			if (!rootNodeFound)
				{
				rootNodeFound = true;
				check(node.motherNode_.ni_ == Def::invalidIdx, "root node has no mother node");
				// 9) root node is all edges
				check(node.isEdgeNode_ == 0x0f0f,"root node is all edges");
				}
			else
				check(node.motherNode_.ni_ != Def::invalidIdx, "non-root node has a mother node");

			if (!node.IsLeaf())
				{
				const TNode& lowChildNode = kdt.Node(node.lowChildNode_);
				const TNode& hiChildNode = kdt.Node(node.hiChildNode_);
				//10) test correct splitting of bounding box into child nodes
				for (size_t i = 0; i < Def::dim; ++i)
					{
					if (i == node.splitDim_)
						{
						check(lowChildNode.corner1_[i] == hiChildNode.corner0_[i]
							&& lowChildNode.corner0_[i] == node.corner0_[i]
							&& hiChildNode.corner1_[i] == node.corner1_[i],
							"correct splitting of bounding box into child nodes at splitDim_");
						}
					else
						{
						check(lowChildNode.corner0_[i] == node. corner0_[i]
							&& hiChildNode.corner0_[i] == node.corner0_[i]
							&& lowChildNode.corner1_[i] == node.corner1_[i]
							&& hiChildNode.corner1_[i] == node.corner1_[i],
							"no splitting of bounding box into child nodes at non-splitDim_");
						}
					}
				// 11) test correct edge node treatment
				std::uint16_t loFlag = ~(1 << (8 + node.splitDim_));
				std::uint16_t hiFlag = ~(1 << node.splitDim_);
				check((node.isEdgeNode_& loFlag) == lowChildNode.isEdgeNode_
					&& (node.isEdgeNode_ & hiFlag) == hiChildNode.isEdgeNode_,
					"correct edge node treatment");
				}
			// 12) all points are within bounding box
			for (Def::TIdxIdx i = node.beginpt_; i.ii_ != node.endpt_.ii_; ++i)
				{
				Def::TPointIdx pi = kdt.Index(i);
				Def::TKDPoint p = kdt.Point(pi);
				check(IsInBox(node.Box(), p), "all points are within bounding box of node");
				}
			} // for all nodes
		
		// 13) for i = 0..m-1 
		//			test if nodeIdx_ entry i really points to the leaf node which contains the single point index i
		for (Def::TIdx i = 0; i < m; ++i)
			{
			Def::TPointIdx pi{ i };

			Def::TNodeIdx ni = kdt.NodeIndex(pi);
			const TNode& node = kdt.Node(ni);
			check(node.IsLeaf(), "nodeIdx_ entries must point to leaf nodes");
			Def::TIdxIdx ii = node.beginpt_;
			check(pi.pi_ == kdt.Index(ii).pi_, "nodeIdx_ entry i points to the leaf node which contains the single point index i");
			}


		Def::TNodeIdxArray testLocateIndices = kdt.LocateOverlappingLeafNodes(Def::TBox{ {-2,-2,-1,-1},{2,2,1,1} });
		Def::TBox box{ {-0.5f,-0.5f,-0.5f,-0.5f},{0.5f,0.5f,0.5f,0.5f} };
		Def::TNodeIdxArray testLocateIndices2 = kdt.LocateOverlappingLeafNodes(box);
		Def::TNodeIdxArray testLocateIndices3 = kdt.LocateOverlappingLeafNodes(Def::TBox{ {-0.005f,-0.005f,-0.005f,-0.005f},{0.005f,0.005f,0.005f,0.005f} });
		std::vector<bool> nodeFound(kdt.Nodes().size(), false);
		for (auto ni : testLocateIndices2)
			{
			nodeFound[ni.ni_] = true;
			check(BoxesOverlap(box, kdt.Node(ni).Box()), "LocateOverlappingLeafNodes: found boxes do overlap");
			}
		for (Def::TIdx i = 0; i < nodeFound.size(); ++i)
			{
			if (!nodeFound[i])
				{
				TNode node = kdt.Node(Def::TNodeIdx{ i });
				if (node.IsLeaf())
					check(!BoxesOverlap(box, node.Box()), "LocateOverlappingLeafNodes: other boxes don't overlap");
				}
			}



		for (TKDTree::TIdx i = 0; i < kdt.Points().size(); ++i)
			{
			std::cout << i << ' ';
			std::cout.flush();
			Def::TKDPoint pt = pts[i];
			std::cout << i << "(" << pt[0] << ',' << pt[1] << ',' << pt[2] << ',' << pt[3] << ")\n";
			Def::TNodeIdx idx = kdt.Locate(TKDTree::TPointIdx{ i });
			TNode nd = kdt.Nodes()[idx.ni_];
			assert(nd.IsLeaf());
			assert(i == kdt.Index()[nd.beginpt_.ii_].pi_);

			assert(idx.ni_ == kdt.Locate(pt).ni_);
			}

		for (const auto& node: kdt.Nodes())
			{
			for (TKDTree::TIdxIdx i = node.beginpt_; i.ii_ != node.endpt_.ii_; ++(i.ii_))
				{
				
				}

			}

		}

	TKDTree::TNearestNeighbors::TNearestNeighbors(TIdx n)
		: i_points_(n, TPointIdx{ invalidIdx }),
		i_nodes_(n, TNodeIdx{ invalidIdx }),
		distances_(n, std::numeric_limits<TReal>::max()) {
		}

	Def::TReal TNode::Distance(const Def::TKDPoint& pt) const
		{
		auto sqr = [](Def::TReal r) {return r * r; };
		Def::TReal rv = 0;
		for (int i = 0; i < pt.size(); ++i)
			{
			if (pt[i] < corner0_[i]) rv += sqr(corner0_[i] - pt[i]);
			if (pt[i] > corner1_[i]) rv += sqr(corner1_[i] - pt[i]);
			}
		return sqrt(rv);
		}
	Def::TReal TNode::Volume() const
		{
		Def::TReal rv = corner1_[0] - corner0_[0];
		for (int i = 1; i < Def::dim; ++i)
			rv *= corner1_[i] - corner0_[i];
		return rv;
		}

	// partition node into nPoints sub-blocks which all have same volume.
	// precondition: point must be within this node bounding box
	// returns list of points within sub-blocks. These points are all at center, except for the sub-block which contains pt
	// rv has nPoints-1 center points and pt
	std::vector<Def::TKDPoint> TNode::Partition(Def::TIdx nPoints, const Def::TKDPoint& pt) const
		{
		std::vector<Def::TKDPoint> rv;
		if (nPoints == 0)
			return rv;
		if (!IsInBox({ corner0_,corner1_ }, pt))
			throw std::runtime_error("TNode::Partition: pt not in bounding box");
		if (nPoints == 1)
			{
			rv.push_back(pt);
			return rv;
			}
		struct S { Def::TIdx has_n; size_t slice_dim;  Def::TKDPoint c0; Def::TKDPoint c1; };
		std::stack<S> work;
		size_t slicedim = 0;
		work.push(S{ nPoints, slicedim, corner0_, corner1_ });
		bool ptTakenCareOf = false;
		while (!(work.empty()))
			{
			S todo = work.top();
			work.pop();
			slicedim = todo.slice_dim;
			if (todo.has_n == 1)
				{
				if (!ptTakenCareOf && IsInBox({ todo.c0, todo.c1 }, pt))
					{
					rv.push_back(pt);
					ptTakenCareOf = true;
					}
				else
					rv.push_back(MidPoint(todo.c0, todo.c1));
				}
			else
				{
				using R = Def::TReal;
				Def::TIdx nlo = todo.has_n / 2;
				Def::TIdx nhi = todo.has_n - nlo;
				R fac = static_cast<R>(nlo) / static_cast<R> (todo.has_n);
				Def::TKDPoint c1lo = todo.c1;
				c1lo[todo.slice_dim] = todo.c0[todo.slice_dim] * (1 - fac) + todo.c1[todo.slice_dim] * fac;
				Def::TKDPoint c0hi = todo.c0;
				c0hi[todo.slice_dim] = c1lo[todo.slice_dim];
				S lo{ nlo, (slicedim + 1) % Def::dim, todo.c0, c1lo };
				S hi{ nhi, lo.slice_dim, c0hi, todo.c1 };
				work.push(lo);
				work.push(hi);
				}
			}
		return rv;
		}


	// #pragma optimize( "", off )
	std::vector<Def::TKDPoint> TNode::RandomPartition(Def::TIdx nPoints, const Def::TKDPoint& pt, const std::function<Def::TReal()>& ranGen) const
		{
		std::vector<Def::TKDPoint> rv;
		rv.reserve(nPoints);
		if (nPoints == 0)
			return rv;
		if (!IsInBox({ corner0_,corner1_ }, pt))
			throw std::runtime_error("TNode::Partition: pt not in bounding box");
		if (nPoints == 1)
			{
			rv.push_back(pt);
			return rv;
			}
		for (size_t i = 0; i < nPoints; ++i)
			{
			Def::TKDPoint thispt;
			for (Def::TIdx j = 0; j < Def::dim; ++j)
				{
				Def::TReal random = ranGen();
				thispt[j] = corner0_[j] * (1 - random) + corner1_[j] * random;
				}
			rv.push_back(thispt);
			}
		return rv;
		}
	} // end namespace