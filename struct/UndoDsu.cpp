struct DSU {
	std::vector<int> parent, siz;
	std::vector<std::array<int, 5> > stk;
	DSU(int n) : parent(n + 1), siz(n + 1, 1) {
		std::iota(parent.begin(), parent.end(), 0);
	}	

	int leader(int x) {
		while (x != parent[x]) {
			x = parent[x];
		}
		return x;
	}

	bool merge(int x, int y, int t) {
		x = leader(x), y = leader(y);
		if (x == y) return false;
		if (siz[x] < siz[y]) {
			std::swap(x, y);
		}

		stk.push_back({t, x, siz[x], y, siz[y]});
		siz[x] += siz[y];
		parent[y] = x;
		return true;
	}

	void undo(int t) {
		while (stk.size() && stk.back()[0] > t) {
			auto &[_, x, sx, y, sy] = stk.back();
			siz[x] = sx;
			parent[x] = x;
			siz[y] = sy;
			parent[y] = y;
			stk.pop_back();
		}
	}

};