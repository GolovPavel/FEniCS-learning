# fenics-learning
The repository for learning fenics-project library

## Docker settings
Start docker in Linux:
`sudo service docker start`

Create new docker container
`docker run -ti -v $(pwd):/home/fenics/shared --name CONTAINER-NAME quay.io/fenicsproject/stable`

The container can be be stopped and started:
`docker stop fenics-container
docker start fenics-container
docker exec -ti -u fenics fenics-container /bin/bash -l`

To see the name and other information of every container you have ever created:
`docker ps -a`
