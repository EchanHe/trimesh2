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
		//if (optimized == true) {
		//	float sum = std::accumulate(attri.begin(), attri.end(), 0.0f);
		//	float mean = sum / nv;
		//	std::vector<float> temp(nv,0);
		//	for (int i = 0; i < nv; i++)
		//		temp[i] = sqr(attri[i] - mean);
		//	float stDev =  sqrt(std::accumulate(temp.begin(), temp.end(), 0.0f)
		//		/ nv);
		//	max = (mean + 2 * stDev > *std::max_element(attri.begin(), attri.end())) ? *std::max_element(attri.begin(), attri.end()) : mean + 2 * stDev;
		//	min = (mean - 2 * stDev < *std::min_element(attri.begin(), attri.end())) ? *std::min_element(attri.begin(), attri.end()) : mean - 2 * stDev;
		//	if (mean - 2 * stDev < *std::min_element(attri.begin(), attri.end()))
		//		min = *std::min_element(attri.begin(), attri.end());
		//	else
		//		min = mean - 2 * stDev;
		//	//max = mean + 3 * stDev;
		//	//min = mean - 3 * stDev;
		//}
		if (optimized == true) {
			std::vector<float> temp;
			temp = attri;
			std::sort(temp.begin(), temp.end());

			float percent = 0.90;
			int portion = (int)(percent*nv); 
			std::vector<float> temp2;
			for (int i = nv - portion; i < portion; i++) {
				temp2.push_back(temp[i]);
			}
			min = *std::min_element(temp2.begin(), temp2.end());
			max = *std::max_element(temp2.begin(), temp2.end());
			//	std::cout << "color max min: " << min << " " << max << std::endl;
		}
		else {
			//auto min = std::min_element(attri.begin(), attri.end());
			//auto ma = std::max_element(attri.begin(), attri.end());
			min = *std::min_element(attri.begin(), attri.end());
			max = *std::max_element(attri.begin(), attri.end());
			//std::cout << *ma;
		}
		float scale = 1 / (max - min);

		for (int i = 0; i < attri.size(); i++) {
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

	void TriMesh::color_shape_index() {
		int nv = vertices.size();
		if (shape_index.size() == 0)
			return;
		if (colors.size() == vertices.size())
			return;
		colors.resize(nv);

		for (int i = 0; i < nv; i++) {
			if(shape_index[i]>0.9)
				colors[i] = Color(0.0, 0.0, 1.0);
			else if (shape_index[i]<-0.9)
				colors[i] = Color(1.0, 0.0, 0.0);
			else 
				colors[i] = Color(0.0, 1.0, 0.0);
		}
	}

	void TriMesh::color_vertex(::std::vector<int> index) {
		int nv = vertices.size();
		if (index.size()==0)
			return;
		if (colors.size() == vertices.size())
			return;
		colors.resize(nv);

		for (int i = 0; i < nv; i++) {
			colors[i] = Color(0.2, 0.2, 0.2);
		}

		for (int i = 0; i < index.size(); i++) {
			colors[index[i]] = Color(1.0, 0.1, 0.1);
		}

	}
}

