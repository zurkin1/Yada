## Start from this Docker image
FROM python:3

## Install in Docker image
RUN pip install numpy
RUN pip install pandas
RUN pip install Cython
RUN pip install tslearn
RUN pip install pymc3

ADD yada /usr/local/yada

#RUN chmod a+x /usr/local/bin/run_model.R

#CMD [ "python", "./my_script.py" ]
## Make Docker container executable
ENTRYPOINT ["python", "/usr/local/yada/run.py"]