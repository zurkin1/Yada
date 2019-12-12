## Start from this Docker image
##Build using:
##docker build -t docker.synapse.org/syn19061090/my-model:version1 .
## --no-cache
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
RUN pip install Cython
RUN pip install tslearn
#RUN pip install pymc3

COPY yada/methods.py /methods.py
COPY yada/diffexp.py /diffexp.py
COPY yada/run.py /run.py
COPY yada/data/Challenge/LM8.csv /LM8.csv
COPY yada/data/Challenge/LM14.csv /LM14.csv
COPY yada/data/Challenge/pure_orig8v2_8.csv /pure_orig8v2_8.csv
COPY yada/data/Challenge/pure_orig8v2_14.csv /pure_orig8v2_14.csv
COPY yada/data/Challenge/ImmunoState.csv /ImmunoState.csv

#RUN mkdir /input/
#RUN chmod a+x /usr/local/bin/run_model.R
#CMD [ "python", "./my_script.py" ]
## Make Docker container executable
ENTRYPOINT ["python", "/r3s4.py"]