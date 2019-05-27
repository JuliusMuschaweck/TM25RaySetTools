#include"KDTree.h"
#include<numeric>
#include<algorithm>
#include<limits>
#include<assert.h>
#include<random>
#include<fstream>
#include<chrono>

#include <iostream>

namespace KDTree
	{
	TKDTree::TKDTree(const TKDPoints & pts)
		: pts_(pts)
		{
		Init();
		}

	TKDTree::TKDTree(TKDPoints && pts)
		: pts_(std::move(pts))
		{
		Init();
		}

	void TKDTree::CreateTree()
		{
		Def::TReal big = std::numeric_limits<Def::TReal>::max();
		// construct root node
		nodes_.reserve(2 * pts_.size() - 1);
		TNode root = RootNode();
		TIdx iRoot = AddNode(root);
		PushWork(iRoot);
		while (!work_.empty())
			SplitNext();
		TIdx ii = 0;
		for (auto i : idx_)
			reverseIdx_[i] = ii++;
		}

	void TKDTree::ShrinkEdgeNodes(TReal fac)
		{
		Def::TReal avgPtsPerDim = pow(static_cast<Def::TReal>(pts_.size()), static_cast<Def::TReal>(1.0 / Def::dim));
		TKDPoint d;
		for (int i = 0; i < dim; ++i)
			d[i] = fac * (boundingBox_.second[i] - boundingBox_.first[i]) / avgPtsPerDim;
		for (auto &nd : nodes_)
			{
			if (nd.IsLeaf() && (nd.isEdgeNode_ != 0x0000))
				{
				TKDPoint pt = pts_[idx_[nd.beginpt_]];
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

	const Def::TKDPoints & TKDTree::Points() const
		{
		return pts_;
		}

	const TKDTree::TNodeArray & TKDTree::Nodes() const
		{
		return nodes_;
		// TODO: insert return statement here
		}

	const Def::TIdxArray & TKDTree::Index() const
		{
		return idx_;
		}

	const Def::TIdxArray & TKDTree::ReverseIndex() const
		{
		return reverseIdx_;
		}

	std::pair<Def::TKDPoint, Def::TKDPoint> TKDTree::BoundingBox() const
		{
		return boundingBox_;
		}

	Def::TIdx TKDTree::Locate(const TKDPoint & pt) const
		{
		const TNode* nd = &(nodes_[0]);
		for (int i = 0; i < dim; ++i)
			{
			if (pt[i] < nd->corner0_[i] || pt[i] > nd->corner1_[i])
				{
				return invalidIdx;
				}
			}
		TIdx idx = 0;
		while (!(nd->IsLeaf()))
			{
			size_t sd = nd->splitDim_;
			if (pt[sd] < nodes_[nd->lowChildNode_].corner1_[sd])
				idx = nd->lowChildNode_;
			else
				idx = nd->hiChildNode_;
			nd = &(nodes_[idx]);
			}
		return idx;
		}

	Def::TIdx TKDTree::Locate(TIdx i) const
		{
		if (i >= pts_.size())
			return invalidIdx;
		TIdx ri = reverseIdx_[i];
		TIdx nb = 0;
		while (!(nodes_[nb].IsLeaf()))
			{
			TIdx dlo = nodes_[nb].lowChildNode_;
			if (ri < nodes_[dlo].endpt_)
				nb = dlo;
			else
				nb = nodes_[nb].hiChildNode_;
			}
		return nb;
		}

	// return the index of the nearest node with >= n points if pt is outside bounding box of all nodes
	Def::TIdx ClosestNodeIndex(const TKDTree::TNodeArray& nodes, const Def::TKDPoint & pt, Def::TIdx n)
		{
		Def::TIdx inode = 0;
		for (;;)
			{
			const TNode& thisnode = nodes[inode];
			const TNode& lonode = nodes[thisnode.lowChildNode_];
			const TNode& hinode = nodes[thisnode.hiChildNode_];
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
			if (jTest < (n-1) && nn.distances_[jTest] < nn.distances_[jTest + 1])
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

	TKDTree::TNearestNeighbors TKDTree::NearestNeighbors(const TKDPoint & pt, TIdx n) const
		{
		if (n >= pts_.size())
			throw std::runtime_error("TKDTree::NearestNeighbors: too many points requested");
		TIdx inode = Locate(pt);
		if (inode == invalidIdx)
			{ // pt is outside the bounding box
			// look for closest box with >= n points
			inode = ClosestNodeIndex(nodes_, pt, n);
			}
		return NearestNeighborsOfNode(inode, pt, n);
		}

	TKDTree::TNearestNeighbors TKDTree::NearestNeighborsOfPoint(TIdx ipoint, TIdx n) const
		{
		if (ipoint >= pts_.size())
			throw std::runtime_error("TKDTree::NearestNeighborsOfPoint: ipoint > npoints");
		TIdx inode = ReverseIndex()[ipoint];
		return NearestNeighborsOfNode(inode, pts_[ipoint], n);
		}

	TKDTree::TNearestNeighbors TKDTree::NearestNeighborsOfNode(TIdx inode, TIdx n) const
		{
		// use center point of node
		const TNode& node = nodes_[inode];
		const TKDPoint& lo = node.corner0_;
		const TKDPoint& hi = node.corner1_;
		TKDPoint pt{ 0.5*(lo[0] + hi[0]),0.5*(lo[1] + hi[1]), 0.5*(lo[2] + hi[2]), 0.5*(lo[3] + hi[3]) };
		return NearestNeighborsOfNode(inode, pt, n);
		}

	TKDTree::TNearestNeighbors TKDTree::NearestNeighborsOfNode(TIdx inode, const TKDPoint& pt, TIdx n) const
		{
		if (inode == invalidIdx)
			throw std::runtime_error("TKDTree::NearestNeighborsOfNode: invalidIdx");
		if (inode >= nodes_.size())
			throw std::runtime_error("TKDTree::NearestNeighborsOfNode: inode out of range");
		// find smallest mother box with n points
		while (nodes_[inode].NPoints() < n)
			inode = nodes_[inode].motherNode_;
		// now inode points to nearest box with >= n points
		// select n nearest within this, as heap
		const TNode& nd = nodes_[inode];
		TNearestNeighbors rv(n); // filled with BIG and invalid indices
		for (TIdx ipt = nd.beginpt_; ipt < nd.endpt_; ++ipt)
			{
			const TKDPoint& thispt = pts_[idx_[ipt]];
			TReal d = Distance(thispt, pt);
			if (d < rv.distances_[0])
				{
				// replace top of heap with new values
				rv.distances_[0] = d;
				rv.i_points_[0] = idx_[ipt];
				rv.i_nodes_[0] = nodeIdx_[idx_[ipt]];
				if (n > 1)
					Heapify(rv, n); // maintain heap structure
				}
			}
		// rv is filled with heap of n candidate points, all from same box
		// now traverse tree, considering only potentially better boxes
		std::stack<TIdx> task;
		task.push(0); // start working with top node
		while (!task.empty())
			{
			TIdx itodo = task.top();
			task.pop();
			if (itodo == inode) // this whole node is taken care of already
				continue;
			const TNode& nTodo = nodes_[itodo];
			if (nTodo.Distance(pt) < rv.distances_[0]) // skip whole node if too far away
				{ // candidate node
				if (nTodo.IsLeaf())
					{ // consider point
					TIdx ipt = idx_[nTodo.beginpt_];
					TReal d = Distance(pts_[ipt], pt);
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

	Def::TReal TKDTree::TotalVolume(const TIdxArray & inodes) const
		{
		TReal rv = 0;
		for (auto i : inodes)
			{
			rv += nodes_[i].Volume();
			}
		return rv;
		}

	void TKDTree::Init()
		{
		idx_.resize(pts_.size());
		std::iota(idx_.begin(), idx_.end(), 0);
		reverseIdx_ = idx_;
		nodeIdx_.resize(pts_.size());
		std::fill(nodeIdx_.begin(), nodeIdx_.end(), invalidIdx);
		for (size_t d = 0; d < Def::dim; ++d)
			ptCoords_[d].reserve(pts_.size());
		Def::TReal big = std::numeric_limits<Def::TReal>::max();
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
		TNode root;
		root.beginpt_ = 0;
		root.endpt_ = static_cast<TIdx>(pts_.size());
		std::tie(root.corner0_, root.corner1_) = KDTree::BoundingBox(pts_);
		Def::TReal avgPtsPerDim = pow(static_cast<Def::TReal>(pts_.size()), static_cast<Def::TReal>(1.0 / Def::dim));
		for (int i = 0; i < Def::dim; ++i)
			{
			Def::TReal d = (root.corner1_[i] - root.corner0_[i]) / avgPtsPerDim * static_cast<Def::TReal>(0.1);
			root.corner0_[i] -= d;
			root.corner1_[i] += d;
			}
		root.motherNode_ = Def::invalidIdx; // root node is the only node w/o mother
		root.lowChildNode_ = Def::invalidIdx; // will be created later
		root.hiChildNode_ = Def::invalidIdx; // will be created later
		root.splitDim_ = 0; // we start with splitting along dimension 1
		std::array<uint16_t, 9> edgeFlags = { 0x0000, 0x0101, 0x0303, 0x0707, 0x0f0f, 0x1f1f, 0x3f3f, 0x7f7f, 0xffff };
		root.isEdgeNode_ = edgeFlags[dim];
		return root;
		}

	Def::TIdx TKDTree::AddNode(const TNode & n)
		{
		// std::lock_guard<std::mutex> lock(mutex_);
		// this is where thread safe node adding may come in, e.g. by locking a mutex
		nodes_.push_back(n);
		return static_cast<TIdx>(nodes_.size()) - 1;
		}

	void TKDTree::SetChildNodes(TIdx mother, TIdx loChild, TIdx hiChild)
		{
		nodes_[mother].lowChildNode_ = loChild;
		nodes_[mother].hiChildNode_ = hiChild;
		}

	void TKDTree::PushWork(TIdx iNode)
		{
		work_.push(iNode);
		}

	Def::TIdx TKDTree::PopWork()
		{
		Def::TIdx rv{ work_.top() };
		work_.pop();
		return rv;
		}

	void TKDTree::SplitNext()
		{
		// retrieve node to work on
		TIdx next = PopWork();
		TNode node = nodes_[next];
		// the coordinates according to which the split shall occur
		#ifndef NDEBUG
		// std::cout << "working on node " << next << '\n';
		#endif
		const std::vector<Def::TReal>& coords = ptCoords_[node.splitDim_];
		// the range of indices to be split
		TIdx nPts = node.endpt_ - node.beginpt_;
		assert(nPts >= 2);
		auto ibegin = idx_.begin() + node.beginpt_;
		auto iend = idx_.begin() + node.endpt_;
		auto imid = ibegin + nPts/2;
		// the comparison function for nth_element
		auto comp = [&coords](TIdx i1, TIdx i2) -> bool
			{
			// return coords[this->idx_[i1]] < coords[this->idx_[i2]];
			return coords[i1] < coords[i2];
			};
		// the actual splitting: the range [ibegin, iend) of idx_ is rearranged
		std::nth_element(ibegin, imid, iend, comp);
		// now  [ibegin, imid) has low values of coords, and [imid, iend) high values.
		// imid points to the split point. But we want to split the bounding box between
		// the split point and its nearest lower neighbor
 		auto closestLoIdx = std::max_element(ibegin, imid, comp);
		TReal closestLoCoord = coords[*closestLoIdx];
		TReal splitPointCoord = coords[*imid];
		TReal BBSplitCoord = (closestLoCoord + splitPointCoord) / 2;
		TNode dlo;
		TNode dhi;
		dlo.beginpt_ = node.beginpt_;
		dlo.endpt_ = node.beginpt_ + nPts / 2;
		dhi.beginpt_ = dlo.endpt_;
		dhi.endpt_ = node.endpt_;
		dlo.corner0_ = dhi.corner0_ = node.corner0_;
		dlo.corner1_ = dhi.corner1_ = node.corner1_;
		dlo.corner1_[node.splitDim_] = dhi.corner0_[node.splitDim_] = BBSplitCoord;
		dlo.motherNode_ = dhi.motherNode_ = next;
		dlo.lowChildNode_ = dhi.lowChildNode_ = dlo.hiChildNode_ = dhi.hiChildNode_ = Def::invalidIdx;
		dlo.splitDim_ = dhi.splitDim_ = (node.splitDim_ + 1) % dim;
		dlo.isEdgeNode_ = dhi.isEdgeNode_ = node.isEdgeNode_;
		// adjust bit field of dim
		std::uint16_t loFlag = 1;
		loFlag = loFlag << (8 + node.splitDim_);
		loFlag = ~loFlag;
		dlo.isEdgeNode_ &= loFlag;
		std::uint16_t hiFlag = 1;
		hiFlag = hiFlag << node.splitDim_;
		hiFlag = ~hiFlag;
		dhi.isEdgeNode_ &= hiFlag;
		TIdx ilo = AddNode(dlo);
		TIdx ihi = AddNode(dhi);
		SetChildNodes(next, ilo, ihi);
		if (!dlo.IsLeaf())
			PushWork(ilo);
		else
			nodeIdx_[dlo.beginpt_] = ilo;
		if (!dhi.IsLeaf())
			PushWork(ihi);	
		else
			nodeIdx_[dhi.beginpt_] = ihi;
		}

	Def::TReal Distance(const Def::TKDPoint & p1, const Def::TKDPoint & p2)
		{
		auto sqr = [](Def::TReal r) {return r * r; };
		Def::TReal rv = 0;
		for (int i = 0; i < p1.size(); ++i)
			{
			rv += sqr(p2[i] - p1[i]);
			}
		return sqrt(rv);
		}

	std::pair<Def::TKDPoint, Def::TKDPoint> BoundingBox(const Def::TKDPoints & pts)
		{
		Def::TReal big = std::numeric_limits<Def::TReal>::max();
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
		return std::make_pair(p0,p1);
		}

	void TestKDTree2D(std::string fn)
		{
		assert(Def::dim == 2);
		size_t npts = 1000*1000;
		Def::TKDPoints pts;
		std::mt19937 gen;
		gen.seed();
		std::uniform_real_distribution<float> dis(-1.0f, 1.0f);
		for (size_t i = 0; i < npts; ++i)
			{
			auto pt = Def::TKDPoint{ 2 * dis(gen), dis(gen) };
			if ((pt[0] * pt[0])/4 + pt[1] * pt[1] < 1)
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
		
		return; 
		std::ofstream f(fn);
		f << "% output of TestKDTree2D\n";
		f << "clear; figure(1); clf; hold on;\n";
		for (auto n : nd)
			{
			if (n.IsLeaf())
				{
				Def::TKDPoint pt = pts[n.beginpt_];
				Def::TKDPoint bblo = n.corner0_;
				Def::TKDPoint bbhi = n.corner1_;
				f << "plot([" << bblo[0] << ',' << bbhi[0] << ',' << bbhi[0] << ',' << bblo[0] << ',' << bblo[0] << "], ";
				f << '[' << bblo[1] << ',' << bblo[1] << ',' << bbhi[1] << ',' << bbhi[1] << ',' << bblo[1] << "],'k');\n";
				f << "scatter(" << pt[0] << "," << pt[1] << ",'kx');\n";
				}
			}
		for (auto n : nd)
			{
			if (n.IsLeaf())
				{
				Def::TKDPoint pt = pts[n.beginpt_];
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
			auto pt1 = pts[nn00.i_points_[i]];
			auto pt2 = pts[nn22.i_points_[i]];
			f << "scatter([" << pt1[0] << ',' << pt2[0] << "],[" << pt1[1] << ',' << pt2[1] << "],'m');\n";
			}
		f << "axis equal;\n";
		}


	void TestKDTree4D()
		{
		assert(Def::dim == 4);
		size_t npts = 1000;
		Def::TKDPoints pts;
		std::mt19937 gen;
		gen.seed();
		std::uniform_real_distribution<float> dis(-1.0f, 1.0f);
		for (size_t i = 0; i < npts; ++i)
			{
			Def::TKDPoint pt;
			std::generate(pt.begin(), pt.end(), [&]() {return dis(gen); });
			if (pt[0] * pt[0] + pt[1] * pt[1] < 1)
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

		for (int i = 0; i < kdt.Points().size(); ++i)
			{
			std::cout << i << ' ';
			std::cout.flush();
			Def::TKDPoint pt = pts[i];
			Def::TIdx idx = kdt.Locate(i);
			TNode nd = kdt.Nodes()[idx];
			assert(nd.IsLeaf());
			assert(i == kdt.Index()[nd.beginpt_]);
			
			assert(idx == kdt.Locate(pt));

			}

		}

	TKDTree::TNearestNeighbors::TNearestNeighbors(TIdx n)
		: i_points_(n, invalidIdx),
		  i_nodes_(n, invalidIdx),
		  distances_(n,std::numeric_limits<TReal>::max())		{
		}

	Def::TReal TNode::Distance(const Def::TKDPoint & pt) const
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
} // end namespace