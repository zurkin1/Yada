## Start from this Docker image
##Build using:
##docker build -t docker.synapse.org/syn19061090/my-model:version1 .
##Push using:
##docker login docker.synapse.org
##docker images
##docker push docker.synapse.org/syn19061090/my-model:version1
##docker run docker.synapse.org/syn19061090/my-model:version1

FROM python:3

## Install in Docker image
RUN pip install numpy
RUN pip install pandas
RUN pip install sklearn
#RUN pip install Cython
#RUN pip install tslearn
#RUN pip install pymc3

COPY yada/run_challenge.py /run_challenge.py
COPY yada/diffexp.py /diffexp.py
COPY yada/data/Challenge/pure.csv /pure.csv
#Files for testing locally.
RUN mkdir /input/
#COPY yada/data/Challenge/input.csv /input/input.csv
#COPY yada/data/EPIC/mix.csv /input/mix.csv
#COPY yada/data/Challenge/example/ds1.csv /input/ds1.csv
#COPY yada/data/Challenge/example/ds2.csv /input/ds2.csv

#RUN chmod a+x /usr/local/bin/run_model.R

#CMD [ "python", "./my_script.py" ]
## Make Docker container executable
ENTRYPOINT ["python", "/run_challenge.py"]