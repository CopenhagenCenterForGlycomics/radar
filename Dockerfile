FROM ubuntu:22.04 AS compile-image

RUN apt-get update
RUN apt-get install -y --no-install-recommends build-essential autoconf python3-setuptools python3-setuptools-whl cython3 python3-wheel  python2-dev python-pip python-setuptools

WORKDIR /radar
ADD . /radar
RUN bash build.sh


FROM ubuntu:22.04 AS runtime-image
RUN apt-get update
RUN apt-get install -y --no-install-recommends python2 python-setuptools

COPY --from=compile-image /usr/local/bin/lfasta /usr/local/bin/lfasta
COPY --from=compile-image /usr/local/bin/radar.py /usr/local/bin/radar.py
COPY --from=compile-image /usr/local/lib/python2.7/dist-packages/Radar* /usr/local/lib/python2.7/dist-packages/

CMD ["/usr/bin/python2", "/usr/local/bin/radar.py"]