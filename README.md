# :closed_book: FEniCS-learning
The repository for learning fenics-project library for university practice.

## Docker settings
Start docker in Linux:
```
sudo service docker start
```
Create new docker container:
```
docker run -ti -v $(pwd):/home/fenics/shared --name CONTAINER-NAME quay.io/fenicsproject/stable
```
The container can be be stopped and started:
```
docker stop fenics-container
docker start fenics-container
docker exec -ti -u fenics fenics-container /bin/bash -l
```
To see the name and other information of every container you have ever created:
```
docker ps -a
```

## Use graphical applications on Linux hosts
This allows X11 applications (e.g. matplotlib plot windows) to be displayed on Linux host systems. To enable this, first run xhost + and then append -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix to the Docker run command. For example, you can run the stable version with:
```
xhost +
docker run -ti -e DISPLAY=$DISPLAY \
   -v /tmp/.X11-unix:/tmp/.X11-unix \
   quay.io/fenicsproject/stable
```
After exiting docker, execute xhost - on the host to restore X settings.

## Running Jupyter notebooks
First of all we run a new Docker container with the jupyter-notebook command specified and the default port 8888 exposed on localhost:
```
docker run --name notebook -w /home/fenics -v $(pwd):/home/fenics/shared -d -p 127.0.0.1:8888:8888 quay.io/fenicsproject/stable 'jupyter-notebook --ip=0.0.0.0'
```
Jupyter creates a unique access token (password) to ensure that only you can access the notebook. To find out this access token run the following command on the host:
```
docker logs CONTAINER-ID
```
Basic two and three-dimensional plotting are available from within the Jupyter notebook.
For matplotlib plotting (2D), open up a new Jupyter notebook, and in the first cell type:
```
%matplotlib inline
```
Execute (Shift-Enter) the cell. In the next cell, we will load in the code from the DOLFIN Python Poisson demo:
```
%load ~/demo/documented/poisson/python/demo_poisson.py
```
Execute (Shift-Enter) the cell. In the same cell, the code from the demo_poisson.py file will be shown. Click in the cell and execute (Shift-Enter) again. A plot of the solution variable u will appear.

For X3DOM plotting (3D), continuing from above, in a new cell type:
```
from IPython.display import HTML
HTML(X3DOM().html(u))
```
Execute (Shift-Enter) the cell. A 3D plot will appear that you can rotate and zoom using the mouse.

## Using ParaView for solution visualisation
ParaView is a powerful tool for visualizing scalar and vector fields, including those computed by FEniCS.
The first step is to export the solution in VTK format:
```
vtkfile = File(’poisson/solution.pvd’)
vtkfile << u
```
The following steps demonstrate how to create a plot of the solution of simple Poisson problem (see samples) in ParaView.
1. Start the ParaView application.
2. Click File–Open... in the top menu and navigate to the directory containing the exported solution. This should be inside a subdirectory named poisson below the directory where the FEniCS Python program was started. Select the file named solution.pvd and then click OK.
3. Click Apply in the Properties pane on the left. This will bring up a plot of the solution.
4. To make a 3D plot of the solution, we will make use of one of ParaView’s many filters. Click Filters–Alphabetical–Warp By Scalar in the top menu and then Apply in the Properties pane on the left. This create an elevated surface with the height determined by the solution value.
5. To show the original plot below the elevated surface, click the little eye icon to the left of solution.pvd in the Pipeline Browser pane on the left. Also click the little 2D button at the top of the plot window to change the visualization to 3D. This lets you interact with the plot by rotating (left mouse button) and zooming (Ctrl + left mouse button).
6. To show the finite element mesh, click on solution.pvd in the Pipeline Browser, navigate to Representation in the Properties pane, and select Surface With Edges. This should make the finite element mesh visible.
7. To change the aspect ratio of the plot, click on WarpByScalar1 in warped plot. We also unclick Orientation Axis Visibility at the bottom of the Properties pane to remove the little 3D axes in the lower left corner of the plot window.
8. Finally, to export the visualization to a file, click File–Save Screen-shot... and select a suitable file name such as poisson.png.

### Animation visualisation in ParaView
To be able to visualize animation solution in an external program such as ParaView, we will save the solution to a file in VTK format in each time step. We do this by first creating a File with the suffix .pvd:
```
vtkfile = File(’some/solution.pvd’)
```
Inside the time loop, we may then append the solution values to this file:
```
vtkfile << (u, t)
```
To visualize the solution, start ParaView, choose File–Open..., open some/solution.pvd, and click Apply in the Properties pane. Click on the play button to display an animation of the solution. To save the animation to a file, click File–Save Animation... and save the file to a desired file format, for example AVI or Ogg/Theora.
