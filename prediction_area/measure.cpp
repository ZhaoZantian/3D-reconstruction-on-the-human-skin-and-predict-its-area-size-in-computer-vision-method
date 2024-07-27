#include "measure.h"

typedef pcl::PointXYZRGB PointT;
typedef pcl::PointCloud<PointT> PointCloud;

double AreaMeasure(string path, double From3DToRealScale)
{
	// fit B-spline surface
//	pcl::on_nurbs::NurbsDataSurface data;
//
//	pcl::visualization::PCLVisualizer viewer("B");
//	viewer.setBackgroundColor(0, 0, 0);
//	viewer.setSize(800, 600);
//	viewer.initCameraParameters();
//	//viewer.addPointCloud<PointT>(cloud, "cloud_cylinder");
//	std::cout << cloud->size() << "points in data set\n" << endl;
//
//	for (unsigned i = 0; i < cloud->size(); i++)
//	{
//		PointT& p = cloud->at(i);
//		data.interior.push_back(Eigen::Vector3d(p.x, p.y, p.z));
//	}
//
//	// parameters
//	unsigned order(3);//B样条曲面的多项式阶
//	unsigned refinement_iterations(3);//求精迭代的次数
//	unsigned mesh_resolution(128);//每个参数方向上的顶点数
//	bool two_dim = true;
//
//	pcl::on_nurbs::FittingSurface::Parameter params;
//	params.interior_smoothness = 0.2;
//	params.interior_weight = 1.0;
//	params.boundary_smoothness = 0.2;
//	params.boundary_weight = 0.0;
//
//	// initialize
//	printf(" surface fitting ...\n");
//	ON_NurbsSurface nurbs = pcl::on_nurbs::FittingSurface::initNurbsPCABoundingBox(order, &data);
//	pcl::on_nurbs::FittingSurface fit(&data, nurbs);
//	// fit.setQuiet (false); // enable/disable debug output
//
//	// mesh for visualization
//	pcl::PolygonMesh mesh;
//	pcl::PointCloud<pcl::PointXYZ>::Ptr mesh_cloud(new pcl::PointCloud<pcl::PointXYZ>);
//	std::vector<pcl::Vertices> mesh_vertices;
//	std::string mesh_id = "mesh_nurbs";
//	pcl::on_nurbs::Triangulation::convertSurface2PolygonMesh(fit.m_nurbs, mesh, mesh_resolution);
//	viewer.addPolygonMesh(mesh, mesh_id);
//	std::cout << "Before refine" << endl;
//
//	for (unsigned i = 0; i < refinement_iterations; i++)
//	{
//		fit.refine(0);
//		if (two_dim)fit.refine(1);
//		fit.assemble(params);
//		fit.solve();
//		pcl::on_nurbs::Triangulation::convertSurface2Vertices(fit.m_nurbs, mesh_cloud, mesh_vertices, mesh_resolution);
//		viewer.updatePolygonMesh<pcl::PointXYZ>(mesh_cloud, mesh_vertices, mesh_id);
//		viewer.spinOnce(1000);
//		std::cout << "refine: " << i << endl;
//	}
//#if 0
//	pcl::PLYWriter writer;
//	writer.write(path, *mesh_cloud, true);
//#else
//	pcl::PCLPointCloud2::Ptr cloud_blob(new pcl::PCLPointCloud2);
//	pcl::toPCLPointCloud2(*mesh_cloud, *cloud_blob);
//	
//	mesh.cloud = *cloud_blob;
//	mesh.polygons = mesh_vertices;
//	pcl::io::savePLYFile(path, mesh);
//#endif
//	while (!viewer.wasStopped())
//	{
//		viewer.spinOnce(5000);
//	}

#if 1
	/* read polydata from mesh*/
	vtkNew<vtkPLYReader> reader;
	const char* path0 = path.c_str();

	reader->SetFileName(path0);
	reader->Update();

	vtkNew <vtkSmoothPolyDataFilter> smooth;
	smooth->SetInputConnection(reader->GetOutputPort());
	smooth->SetNumberOfIterations(20);
	smooth->SetRelaxationFactor(1.0);
	smooth->FeatureEdgeSmoothingOn();
	smooth->BoundarySmoothingOn();
	smooth->Update();
#else
	vtkNew<vtkPolyData> polys;
	pcl::io::pointCloudTovtkPolyData(*mesh_cloud, polys);

	vtkNew<vtkSurfaceReconstructionFilter>surf;
	surf->SetInputData(polys);
	surf->Update();
	vtkNew<vtkContourFilter>contour;
	contour->SetInputConnection(surf->GetOutputPort());
	contour->SetValue(0, 0.0);
	contour->Update();
#endif
	vtkNew<vtkPolyDataMapper> polysmap;
	polysmap->SetInputData(smooth->GetOutput());
	//polysmap->SetInputConnection(contour->GetOutputPort());

	vtkNew<vtkActor> actor;
	actor->SetMapper(polysmap);

	vtkNew<vtkRenderer> render;
	render->AddActor(actor);
	render->SetBackground(0.5, 0.6, 0.7);

	vtkNew<vtkRenderWindow> renderwindow;
	renderwindow->AddRenderer(render);
	renderwindow->SetSize(640, 480);
	renderwindow->SetWindowName("polyData Structure");
	renderwindow->Render();

	vtkNew<vtkRenderWindowInteractor> renderwindowInterac;
	renderwindowInterac->SetRenderWindow(renderwindow);
	renderwindowInterac->Start();

	vtkNew<vtkFillHolesFilter> fillHolesFilter;
	fillHolesFilter->SetInputConnection(smooth->GetOutputPort());
	//fillHolesFilter->SetInputConnection(contour->GetOutputPort());
	fillHolesFilter->Update();

	vtkNew<vtkTriangleFilter> triFilter;
	triFilter->SetInputConnection(fillHolesFilter->GetOutputPort());
	triFilter->Update();

	vtkNew<vtkMassProperties> polygonProperties;
	polygonProperties->SetInputConnection(triFilter->GetOutputPort());
	polygonProperties->Update();

	double area = polygonProperties->GetSurfaceArea() * From3DToRealScale;
	cout << "area = "<< area << endl;
	return area;
}

//distance
//double Scale(pcl::PointCloud<PointT>::Ptr cloud)
//{
//	//segment
//	//pcl::ModelCoefficients::Ptr coefficients(new pcl::ModelCoefficients); //存储输出的模型的系数
//	//pcl::PointIndices::Ptr inliers(new pcl::PointIndices);	//存储内点，使用的点
//	////创建分割对象
//	//pcl::SACSegmentation<PointT> seg;
//	////可选设置
//	//seg.setOptimizeCoefficients(true);
//	////必须设置
//	//seg.setModelType(pcl::SACMODEL_CIRCLE2D); //设置模型类型，检测平面
//	//seg.setMethodType(pcl::SAC_RANSAC);		//设置方法 聚类或随机样本一致性
//	//seg.setDistanceThreshold(0.1);
//	//seg.setInputCloud(cloud);
//	//seg.segment(*inliers, *coefficients);
//
//	//pcl::PointCloud <PointT>::Ptr cloud_circle(new pcl::PointCloud <PointT>);
//	//pcl::copyPointCloud(*cloud, inliers->indices, *cloud_circle);
//
//	//pcl::visualization::PCLVisualizer viewer("Cloud viewer");
//	//pcl::visualization::PointCloudColorHandlerCustom<PointT> cloudcolor_r(cloud_circle, 255, 0, 0);
//
//	//viewer.addCoordinateSystem();
//	//viewer.setBackgroundColor(0, 0, 0);
//	//viewer.addPointCloud(cloud_circle, cloudcolor_r, "cloud_plane");
//
//	//while (!viewer.wasStopped())
//	//{
//	//	viewer.spinOnce(5000);
//	//}
//
//	pcl::search::KdTree<PointT>::Ptr tree(new pcl::search::KdTree<PointT>);
//	tree->setInputCloud(cloud);
//
//	std::vector<pcl::PointIndices> cluster_indices;
//	pcl::EuclideanClusterExtraction<PointT> ec;//创建欧式聚类分割对象
//	ec.setClusterTolerance(0.8); //近邻搜索的搜索半径 /1204 - 0.8; 0221 - 0.02; 0228 - 
//	ec.setMinClusterSize(100); //最小聚类尺寸
//	ec.setMaxClusterSize(100000);
//	ec.setSearchMethod(tree);
//	ec.setInputCloud(cloud);
//	ec.extract(cluster_indices);
//
//	std::vector<Eigen::Vector4f> Eucluextra; //用于储存欧式分割后的点云质心
//	for (std::vector<pcl::PointIndices>::const_iterator it = cluster_indices.begin(); it != cluster_indices.end(); ++it)
//	{
//		pcl::PointCloud<PointT>::Ptr cloud_cluster(new pcl::PointCloud<PointT>);
//		for (std::vector<int>::const_iterator pit = it->indices.begin(); pit != it->indices.end(); pit++)
//			cloud_cluster->points.push_back(cloud->points[*pit]);
//		cloud_cluster->width = cloud_cluster->points.size();
//		cloud_cluster->height = 1;
//		cloud_cluster->is_dense = true;
//		Eigen::Vector4f centroid;  //质心 
//		pcl::compute3DCentroid(*cloud_cluster, centroid); // 计算质心
//		Eucluextra.push_back(centroid);
//	}
//	int n = Eucluextra.size();
//	double result = DBL_MAX;
//	vector<double> min_5;
//	for(int i = 0; i< n - 1; ++i)
//	{
//		for (int j = i + 1; j < n - 1; ++j) 
//		{
//			Eigen::Vector4f dis = Eucluextra[i] - Eucluextra[j];
//			min_5.push_back(sqrt(pow(dis[0], 2) + pow(dis[1], 2) + pow(dis[3], 2)));
//		}
//	}
//	sort(min_5.begin(), min_5.end(), less<double>());
//	double avg_min5 = (min_5[0] + min_5[1] + min_5[2] + min_5[3] + min_5[4]) / 5;
//
//	return RealDis / avg_min5;
//}
	