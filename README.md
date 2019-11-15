Registration Library

This is a registration library using cpu, and gpu.
The aim is to achieve real-time deformable registration.\\


Desing decisions:\\
* Library is header only. This simplify template creation. Template class definition short and only brief descriptions. Functions names intuitive.\\
* ITK core for image interface (reading and writing).\\
* Eigen and ViennaCL as core matrix libraries.\\
* Mesh grid as part of image class. To shared with transformations.

