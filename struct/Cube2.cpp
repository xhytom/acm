//        1  2
//        3  4
// 13 14  5  6 17 18 21 22
// 15 16  7  8 19 20 23 24
//        9 10
//       11 12

struct Cube {
	std::array<int, 25> a;
	Cube (std::array<int, 25> b) {
		a = b;
	}	
	
	void R1234() {
		std::tie(a[1], a[2], a[3], a[4], a[13], a[14], a[5], a[6], a[17], a[18], a[21], a[22]) =
 std::make_tuple(a[3], a[1], a[4], a[2], a[5], a[6], a[17], a[18], a[21], a[22], a[13], a[14]);
	}
	void RIGHT() {
		std::tie(a[1], a[2], a[3], a[4], a[13], a[14], a[5], a[6], a[17], a[18], a[21], a[22]) =
 std::make_tuple(a[3], a[1], a[4], a[2], a[5], a[6], a[17], a[18], a[21], a[22], a[13], a[14]);
		std::tie(a[9], a[10], a[11], a[12], a[15], a[16], a[7], a[8], a[19], a[20], a[23], a[24]) =
 std::make_tuple(a[10], a[12], a[9], a[11], a[7], a[8], a[19], a[20], a[23], a[24], a[15], a[16]);
	}
	void DOWN() {
		std::tie(a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9], a[10], a[11], a[12], a[13], 
				a[14], a[15], a[16], a[17], a[18], a[19], a[20], a[21], a[22], a[23], a[24]) =
 std::make_tuple(a[24], a[23], a[22], a[21], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], 
 				a[15], a[13], a[16], a[14], a[18], a[20], a[17], a[19], a[12], a[11], a[10], a[9]);

	}
	bool ok() {
		return a[1] == a[2] && a[2] == a[3] && a[3] == a[4]
			&& a[5] == a[6] && a[6] == a[7] && a[7] == a[8]
			&& a[9] == a[10] && a[10] == a[11] && a[11] == a[12]
			&& a[13] == a[14] && a[14] == a[15] && a[15] == a[16]
			&& a[17] == a[18] && a[18] == a[19] && a[19] == a[20]
			&& a[21] == a[22] && a[22] == a[23] && a[23] == a[24];
	}
};
