# fundamental introduction
These are my undergraduate graduate project which mainly constructs a method to reconstruct a series of images(especially people's wounds) in 3D level, then using our ways to predict and calculate the area of wounds. We know that the surface wound of the human skin is uneven, and it is difficult to directly calculate the direct contact measurement or traditional algorithm. Therefore, we plan to better predict the area of ​​the wound in combination with computer vision.
# pipeline
The main research problem is the improvement of the three -dimensional reconstruction methods of human surface skin. pass The open source library of the two computer vision of OpenMVG and OpenMVS, after the installation is completed, conducts Windows and Linux experiments through my computer and remote computer room hosts. You need to use the camera you purchased in advance to take the personal body in the laboratory. Because the special structure of the camera can capture pictures of a certain size wound of the human arm, use multiple groups of photographic data sets as input, and resume from the images of multiple perspectives through OpenMVG. The three - dimensional structure of the scene realizes three - dimensional reconstruction and camera positioning functions. Then, through OpenMVS, it is responsible for the three -dimensional matching, in - depth estimation, and three - dimensional reconstruction of the image of multiple perspectives to generate three - dimensional models with geometric and texture information. Finally, the formal file generated by the MeshLab simulation software was reconstructed from densely clouds and completed texture mapping. We plan to adjust the grid and change the sample size, modify the camera parameters to perform experiments, and check the verification and improvement effect.
# tools
useful tools/platforms:Opencv, CMake, [OpenMVG/OpenMVS](https://github.com/openMVG/openMVG/tree/develop?tab=readme-ov-file#openmvg-open-multiple-view-geometry)

# ![MeshLab Logo](src/meshlab/images/eye64.png) MeshLab

[![BuildMeshLab](https://github.com/cnr-isti-vclab/meshlab/actions/workflows/BuildMeshLab.yml/badge.svg)](https://github.com/cnr-isti-vclab/meshlab/actions/workflows/BuildMeshLab.yml)
