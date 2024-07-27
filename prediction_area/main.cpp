#include <iostream>
#include <thread>
#include <vector>

#include "measure.h"
//#include "Select.h"

int main() {
	//1. read pointcloud & segment
	std::string filename = R"(D:\wound\0424_wound\scene_dense_mesh_refine_wound_texture.ply)";
	std::string calibfilename = R"(D:\wound\0424_wound\scene_dense_mesh_refine_seg_texture.ply)";
	
	double pi = 3.1415926;
	double r = 10;//9.7 ½ðÊô»·

	//PointSegment(cloud);

	//double From3DToRealScale = Scale(cloud_calib);

	double area = AreaMeasure(calibfilename, 1);
	double From3DToRealScale = pi * pow(r,2) / area;

	cout << "Scale = " << From3DToRealScale << endl;
	AreaMeasure(filename, From3DToRealScale);

	return (0);
}
