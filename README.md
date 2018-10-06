# fenics-learning
The repository for learning fenics-project library

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
