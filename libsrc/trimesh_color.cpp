#include "TriMesh.h"
#include <algorithm>
#include <numeric>

namespace trimesh {
	void TriMesh::color_vertex(::std::vector<float> attri , bool optimized /*= true*/) {
	


		int nv = vertices.size();
		if (nv != attri.size())
			return;
		if (colors.size() == vertices.size())
			return;
		colors.resize(vertices.size());
		float max, min;
		if (optimized == true) {
			float sum = std::accumulate(attri.begin(), attri.end(), 0.0f);
			float mean = sum / nv;
			std::vector<float> temp(nv,0);
			for (int i = 0; i < nv; i++)
				temp[i] = sqr(attri[i] - mean);
			float stDev =  sqrt(std::accumulate(temp.begin(), temp.end(), 0.0f)
				/ nv);
			max = (mean + 2 * stDev > *std::max_element(attri.begin(), attri.end())) ? *std::max_element(attri.begin(), attri.end()) : mean + 2 * stDev;
			min = (mean - 2 * stDev < *std::min_element(attri.begin(), attri.end())) ? *std::min_element(attri.begin(), attri.end()) : mean - 2 * stDev;
			if (mean - 2 * stDev < *std::min_element(attri.begin(), attri.end()))
				min = *std::min_element(attri.begin(), attri.end());
			else
				min = mean - 2 * stDev;
			//max = mean + 3 * stDev;
			//min = mean - 3 * stDev;
		}
		else {
			//auto min = std::min_element(attri.begin(), attri.end());
			//auto ma = std::max_element(attri.begin(), attri.end());
			min = *std::min_element(attri.begin(), attri.end());
			max = *std::max_element(attri.begin(), attri.end());
			//std::cout << *ma;
		}
		float scale = 1 / (max - min);

		for (int i = 0; i < vertices.size(); i++) {
			if (attri[i] > max)
				attri[i] = max;
			if (attri[i] < min)
				attri[i] = min;


			float intensity = scale * (attri[i] - min);
			if (intensity < 0.5)
				colors[i] = Color(1 - intensity * 2, intensity * 2, 0.0f);
			else
				colors[i] = Color(0.0f, 1 - (intensity - 0.5f) * 2, (intensity - 0.5f) * 2);
		}
	}
}
