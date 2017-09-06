/*
Szymon Rusinkiewicz
Princeton University

TriMesh_stats.cc
Computation of various statistics on the mesh.
*/

#include "TriMesh.h"
#include <algorithm>
#include <numeric>
using namespace std;


namespace trimesh {

// Compute a variety of statistics.  Takes a type of statistic to compute,
// and what to do with it.
float TriMesh::stat(StatOp op, StatVal val)
{
	vector<float> vals;

	switch (val) {
		case STAT_VALENCE: {
			need_neighbors();
			int nv = vertices.size();
			for (int i = 0; i < nv; i++)
				vals.push_back((float) neighbors[i].size());
			break;
		}
		case STAT_FACEAREA: {
			need_faces();
			int nf = faces.size();
			for (int i = 0; i < nf; i++)
				vals.push_back(len(trinorm(i)));
			break;
		}
		case STAT_ANGLE: {
			need_faces();
			int nf = faces.size();
			for (int i = 0; i < nf; i++)
				for (int j = 0; j < 3; j++)
					vals.push_back(cornerangle(i, j));
			break;
		}
		case STAT_DIHEDRAL: {
			need_across_edge();
			int nf = faces.size();
			for (int i = 0; i < nf; i++)
				for (int j = 0; j < 3; j++) {
					if (across_edge[i][j] < 0)
						continue;
					vals.push_back(dihedral(i, j));
				}
			break;
		}
		case STAT_EDGELEN: {
			need_faces();
			int nf = faces.size();
			for (int i = 0; i < nf; i++)
				for (int j = 0; j < 3; j++)
					vals.push_back(dist(
						vertices[faces[i][j]],
						vertices[faces[i][(j+1)%3]]));
			break;
		}
		case STAT_X: {
			int nv = vertices.size();
			for (int i = 0; i < nv; i++)
				vals.push_back(vertices[i][0]);
			break;
		}
		case STAT_Y: {
			int nv = vertices.size();
			for (int i = 0; i < nv; i++)
				vals.push_back(vertices[i][1]);
			break;
		}
		case STAT_Z: {
			int nv = vertices.size();
			for (int i = 0; i < nv; i++)
				vals.push_back(vertices[i][2]);
			break;
		}
		case STAT_MEAN_CURV: {
			int nv = vertices.size();
			for (int i = 0; i < nv; i++)
				vals.push_back(mean_curv[i]);
			break;
		}
		case STAT_GAUS_CURV: {
			int nv = vertices.size();
			for (int i = 0; i < nv; i++)
				vals.push_back(gaus_curv[i]);
			break;
		}
		case STAT_SDF: {
			int nv = vertices.size();
			for (int i = 0; i < nv; i++)
				vals.push_back(sdf[i]);
			break;
		}
		default:
			return 0.0f;
	}

	int n = vals.size();
	if (!n)
		return 0.0f;

	switch (op) {
		case STAT_MIN:
			return *min_element(vals.begin(), vals.end());

		case STAT_MAX:
			return *max_element(vals.begin(), vals.end());

		case STAT_MEANABS:
			for (int i = 0; i < n; i++)
				if (vals[i] < 0.0f)
					vals[i] = -vals[i];
			// Fall through
		case STAT_MEAN:
			return accumulate(vals.begin(), vals.end(), 0.0f) / n;

		case STAT_RMS:
			for (int i = 0; i < n; i++)
				vals[i] *= vals[i];
			return sqrt(accumulate(vals.begin(), vals.end(), 0.0f)
				/ n);

		case STAT_MEDIAN:
			if (n & 1) {
				nth_element(vals.begin(),
					    vals.begin() + n/2,
					    vals.end());
				return vals[n/2];
			} else {
				nth_element(vals.begin(),
					    vals.begin() + n/2 - 1,
					    vals.end());
				float tmp = vals[n/2 - 1];
				nth_element(vals.begin(),
					    vals.begin() + n/2,
					    vals.end());
				return 0.5f * (tmp + vals[n/2]);
			}

		case STAT_STDEV: {
			float mean = accumulate(vals.begin(), vals.end(), 0.0f)
				/ n;
			for (int i = 0; i < n; i++)
				vals[i] = sqr(vals[i] - mean);
			return sqrt(accumulate(vals.begin(), vals.end(), 0.0f)
				/ n);
		}

		case STAT_TOTAL:
			return accumulate(vals.begin(), vals.end(), 0.0f);
	}

	return 0.0f; // Can't happen, I hope.
}


// A characteristic "feature size" for the mesh.  Computed as an approximation
// to the median edge length
float TriMesh::feature_size()
{
	need_faces();
	if (faces.empty())
		return 0.0f;

	int nf = faces.size();
	int nsamp = min(nf / 2, 333);

	vector<float> samples;
	samples.reserve(nsamp * 3);

	for (int i = 0; i < nsamp; i++) {
		// Quick 'n dirty portable random number generator
		static unsigned randq = 0;
		randq = unsigned(1664525) * randq + unsigned(1013904223);

		int ind = randq % nf;
		const point &p0 = vertices[faces[ind][0]];
		const point &p1 = vertices[faces[ind][1]];
		const point &p2 = vertices[faces[ind][2]];
		samples.push_back(dist2(p0,p1));
		samples.push_back(dist2(p1,p2));
		samples.push_back(dist2(p2,p0));
	}
	nth_element(samples.begin(),
		    samples.begin() + samples.size()/2,
		    samples.end());
	return sqrt(samples[samples.size()/2]);
}

void TriMesh::remove_outlier(StatVal val) {
	vector <float>* vals;
	//switch (val) {
	//case STAT_MEAN_CURV: {
	//	for (int i; i < mean_curv.size(); i++) {
	//		if (mean_curv[i] > max)
	//			mean_curv[i] = max;
	//		if (mean_curv[i] < min)
	//			mean_curv[i] = min;
	//	}
	//	break;
	//}
	//default:
	//	return;
	//}

	


	float mean = stat(STAT_MEAN, val);
	float std = stat(STAT_STDEV, val);

	float max = mean + 2 * std;
	float min = mean - 2 * std;

	switch (val) {
		case STAT_MEAN_CURV: {
			for (int i=0; i < mean_curv.size(); i++) {
				if (mean_curv[i] > max)
					mean_curv[i] = max;
				if (mean_curv[i] < min)
					mean_curv[i] = min;
			}
			break;
		}
		case STAT_GAUS_CURV: {
			for (int i=0; i < gaus_curv.size(); i++) {
				if (gaus_curv[i] > max)
					gaus_curv[i] = max;
				if (gaus_curv[i] < min)
					gaus_curv[i] = min;
			}
			break;
		}
		default:
			return;
	}

}

void TriMesh::need_feature_points_curv_sdf()
{
	if (sdf.size() == 0 || mean_curv.size() ==0 || gaus_curv.size() == 0) {
		return ;
	}

	float stdevN = 1;
	int nv = vertices.size();
	//index.clear();
	vector<int> index;

	std::vector<float> temp;
	float percent = 0.95;
	int portion = (int)(percent*nv);

	//calculate mean curvature
	float mean_meanCurv = stat(STAT_MEAN, STAT_MEAN_CURV);
	float std_meanCurv = stat(STAT_STDEV, STAT_MEAN_CURV);
	temp = mean_curv;
	std::sort(temp.begin(), temp.end());
	float mean_curv_lower_thre = temp[1 - portion];
	float mean_curv_upper_thre = temp[portion];
	//calculate gaussian curvature
	float mean_gausCurv = stat(STAT_MEAN, STAT_GAUS_CURV);
	float std_gausCurv = stat(STAT_STDEV, STAT_GAUS_CURV);
	temp = gaus_curv;
	std::sort(temp.begin(), temp.end());
	float gaus_curv_upper_thre = temp[portion];
	float gaus_curv_lower_thre = temp[1-portion];
	//calculate SDF
	//if (sdf.size() != 0) {
	float percent2 = 0.95;
	int portion2 = (int)(percent2*nv);


	float mean_sdf = stat(STAT_MEAN, STAT_SDF);
	float std_sdf = stat(STAT_STDEV, STAT_SDF);
	float median_sdf = stat(STAT_MEDIAN, STAT_SDF);
	float max_sdf = stat(STAT_MAX, STAT_SDF);
	float min_sdf = stat(STAT_MIN, STAT_SDF);
	temp = sdf;
	std::sort(temp.begin(), temp.end());
	float sdf_upper_thre = temp[portion2];
	float sdf_lower_thre = temp[1- portion2];
	//}


	float min_meanCurv = mean_meanCurv - stdevN*std_meanCurv;
	float max_meanCurv = mean_meanCurv + stdevN*std_meanCurv;

	float min_gausCurv = mean_gausCurv - stdevN*std_gausCurv;
	float max_gausCurv = mean_gausCurv + stdevN*std_gausCurv;

	/*float min_sdf = mean_sdf - stdevN*std_sdf;
	float max_sdf = mean_sdf + stdevN*std_sdf;*/

	
	std::vector<vec> results;
	for (int i = 0; i < nv; i++) {

		//bool isPeak = (mean_curv[i] < 0 && mean_curv[i] < min_meanCurv) &&
		//	(gaus_curv[i] > 0 && gaus_curv[i] > max_gausCurv) &&
		//	(sdf[i]>sdf_thre);

		bool isPeak = (mean_curv[i] < 0 && mean_curv[i] < mean_curv_lower_thre) &&
			(gaus_curv[i] > 0 && gaus_curv[i] > gaus_curv_upper_thre);
		//&&(sdf[i]>sdf_upper_thre);

		bool isNarrow = sdf[i] < sdf_lower_thre;
		bool isThick = sdf[i] > sdf_upper_thre;
		bool isRavine = (mean_curv[i] > 0 && mean_curv[i] >mean_curv_upper_thre) &&
			(gaus_curv[i] <=0 && gaus_curv[i]<gaus_curv_lower_thre);

		//bool isBeak = sdf[i]<mean_sdf + 2*std_sdf && sdf[i]> mean_sdf - 2*std_sdf;
		bool isFeather = (sdf[i] < sdf_lower_thre && isPeak) && (sdf[i] > sdf_upper_thre && isRavine);


		//if (isPeak) {
		//	results.push_back(vertices[i]);
		//	index.push_back(i);
		//}
			
		if (isThick && isPeak) {
			results.push_back(vertices[i]);
			index.push_back(i);
		}

	}

	landmarks = results;
	landmarkID = index;

	cout << endl << "Find " << landmarks.size() << " feature points using the " << 100 * percent << " % of the SDF, Mean Curv and Gaussian Curv" << endl;
	return;

}

void TriMesh::need_feature_points_curv_sdf_beak()
{
	if (sdf.size() == 0 || mean_curv.size() == 0 || gaus_curv.size() == 0) {
		return;
	}

	float stdevN = 1;
	int nv = vertices.size();
	//index.clear();
	vector<int> index;

	std::vector<float> temp;
	float percent = 0.90;
	int portion = (int)(percent*nv);

	

	//calculate mean curvature
	float mean_meanCurv = stat(STAT_MEAN, STAT_MEAN_CURV);
	float std_meanCurv = stat(STAT_STDEV, STAT_MEAN_CURV);
	temp = mean_curv;
	std::sort(temp.begin(), temp.end());
	float mean_curv_lower_thre = temp[nv - portion];
	float mean_curv_upper_thre = temp[portion];
	//calculate gaussian curvature
	float mean_gausCurv = stat(STAT_MEAN, STAT_GAUS_CURV);
	float std_gausCurv = stat(STAT_STDEV, STAT_GAUS_CURV);
	temp = gaus_curv;
	std::sort(temp.begin(), temp.end());
	float gaus_curv_upper_thre = temp[portion];
	float gaus_curv_lower_thre = temp[nv - portion];
	//calculate SDF
	//if (sdf.size() != 0) {
	float percent2 = 0.90;
	int portion2 = (int)(percent2*nv);


	float mean_sdf = stat(STAT_MEAN, STAT_SDF);
	float std_sdf = stat(STAT_STDEV, STAT_SDF);
	float median_sdf = stat(STAT_MEDIAN, STAT_SDF);
	/*float max_sdf = stat(STAT_MAX, STAT_SDF);
	float min_sdf = stat(STAT_MIN, STAT_SDF);*/
	temp = sdf;
	std::sort(temp.begin(), temp.end());
	float sdf_upper_thre = temp[portion2];
	float sdf_lower_thre = temp[nv - portion2];
	//}
	cout << mean_gausCurv << "  " << std_gausCurv << endl;
	cout << mean_meanCurv << "  " << std_meanCurv << endl;
	cout << mean_curv_lower_thre << "  " << mean_curv_upper_thre << endl;
	cout << gaus_curv_lower_thre << "  " << gaus_curv_upper_thre << endl;
	cout << sdf_lower_thre << "  " << sdf_upper_thre << endl;


	float min_meanCurv = mean_meanCurv - stdevN*std_meanCurv;
	float max_meanCurv = mean_meanCurv + stdevN*std_meanCurv;

	float min_gausCurv = mean_gausCurv - stdevN*std_gausCurv;
	float max_gausCurv = mean_gausCurv + stdevN*std_gausCurv;

	float min_sdf = mean_sdf - stdevN*std_sdf;
	float max_sdf = mean_sdf + stdevN*std_sdf;


	std::vector<vec> results;
	for (int i = 0; i < nv; i++) {

		bool isPit= (mean_curv[i] < 0 && mean_curv[i]< mean_curv_lower_thre) && (gaus_curv[i]<0);
		
		bool isPeak = (mean_curv[i] > 0);// && (gaus_curv[i] > 0);

		bool isRavine = (mean_curv[i] > 0 && mean_curv[i] >mean_curv_upper_thre) &&
			(gaus_curv[i] < 0 && gaus_curv[i]<gaus_curv_lower_thre);
		//&&(sdf[i]>sdf_upper_thre);

		bool isNarrow = sdf[i] < sdf_lower_thre;
		bool isThick = sdf[i] > sdf_upper_thre;
		bool medium_size = sdf[i] < max_sdf && sdf[i]>min_sdf;



		//bool isBeak = sdf[i]<mean_sdf + 2*std_sdf && sdf[i]> mean_sdf - 2*std_sdf;
		bool isFeather = (sdf[i] < sdf_lower_thre && isPeak) && (sdf[i] > sdf_upper_thre && isRavine);


		//if (isPeak) {
		//	results.push_back(vertices[i]);
		//	index.push_back(i);
		//}

		if (isPit && medium_size) {
			results.push_back(vertices[i]);
			index.push_back(i);
		}
		else if(isPeak && isNarrow) {
			landmarks_2.push_back(vertices[i]);
		}

	}

	landmarks = results;
	landmarkID = index;
	cout<<"Use "<<100 * percent << " % of the SDF, Mean Curv and Gaussian Curv" << endl;
	cout << endl << "Find landmark 1: " << landmarks.size() << endl;
	if (landmarks_2.size() > 0) {
		cout << endl << "Find landmark 2: " << landmarks_2.size() << endl;
	}
	
	return;

}

void TriMesh::writeLandMarks() {
	std::cout << "Writing curvature to _LM.csv" << std::endl;
	std::ofstream file;
	char str[80];
	strcpy(str, filename);
	strcat(str, "_LM.csv");
	file.open(str);
	file << "x , y ,z" << "\n";
	for (int i = 0; i < landmarks.size(); i++) {
		file << landmarks[i][0] << "," << landmarks[i][1] << ","<<landmarks[i][2] << "\n";
	}
	file.close();
}

void TriMesh::need_k_mean(int k) {

	int nLM = landmarks_2.size();
	//int nLM = vertices.size();
	
	vector<int> lm_id(nLM);
	std::iota(lm_id.begin(), lm_id.end(), 0);
	random_shuffle(lm_id.begin(), lm_id.end());
	vector<int> k_id(lm_id.begin(), lm_id.begin()+k);
	vector<vector<int>> cluster_id(k);
	vector<vector<vec>> cluster_vertices(k);
	vector<vec> cluster_mean(k);
	for (int i = 0; i < k; i++) {
		cluster_mean[i] = landmarks_2[k_id[i]];
	}
	float obj = 0;
	//for (int i = 0; i < nLM; i++) {
	//	int final_id = 0;
	//	float min_d;
	//	for (int j = 0; j < k; j++) {
	//		float d=dist(vertices[i], cluster_mean[j]);
	//		if (j == 0) {
	//			final_id = 0;
	//			min_d = d;
	//		}
	//		else if (d < min_d) {
	//			final_id = j;
	//			min_d = d;
	//		}
	//	}
	//	obj += min_d;
	//	cluster_id[final_id].push_back(i);
	//	cluster_vertices[final_id].push_back(vertices[i]);
	//}
	//for (int i = 0; i < k; i++) {
	//	float sumX = 0, sumY = 0, sumZ = 0;
	//	int size_per_cluster = cluster_vertices[i].size();
	//	for (int j = 0; j < size_per_cluster; j++) {
	//		sumX += cluster_vertices[i][j][0];
	//		sumY += cluster_vertices[i][j][1];
	//		sumZ += cluster_vertices[i][j][2];
	//	}
	//	cluster_mean[i] = vec(sumX / size_per_cluster, sumY / size_per_cluster, sumZ / size_per_cluster);

	//}
	float prev_obj = obj;
	do {
		//init
		prev_obj = obj;
		obj = 0;
		for (int i = 0; i < k; i++) {
			cluster_vertices[i].clear();
		}
		//end init
		for (int i = 0; i < nLM; i++) {
			int final_id = 0;
			float min_d;
			for (int j = 0; j < k; j++) {
				//float d = dist(vertices[i], cluster_mean[j]);
				float d = dist(landmarks_2[i], cluster_mean[j]);

				if (j == 0) {
					final_id = 0;
					min_d = d;
				}
				else if (d < min_d) {
					final_id = j;
					min_d = d;
				}
			}
			obj += min_d;
			cluster_id[final_id].push_back(i);
			cluster_vertices[final_id].push_back(landmarks_2[i]);
		}
		for (int i = 0; i < k; i++) {
			float sumX = 0, sumY = 0, sumZ = 0;
			int size_per_cluster = cluster_vertices[i].size();
			for (int j = 0; j < size_per_cluster; j++) {
				sumX += cluster_vertices[i][j][0];
				sumY += cluster_vertices[i][j][1];
				sumZ += cluster_vertices[i][j][2];
			}
			cluster_mean[i] = vec(sumX / size_per_cluster, sumY / size_per_cluster, sumZ / size_per_cluster);
		}
		cout << obj << " " << prev_obj << endl;
	} while (obj != prev_obj);
	clustered_LM.resize(k);
	clustered_LM_mean.resize(k);
	for (int i = 0; i < k; i++) {
		clustered_LM[i] = cluster_vertices[i];
		clustered_LM_mean[i] = cluster_mean[i];
		//for (int j = 0; j < cluster_vertices[i].size(); j++) {
		//	clustered_LM[i][j] = cluster_vertices[i][j];
		//}
	}
}

}; // namespace trimesh
