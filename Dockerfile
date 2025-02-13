FROM ubuntu:latest AS compile-image

RUN apt-get update
RUN apt-get install -y --no-install-recommends build-essential autoconf python3-setuptools python3-setuptools-whl cython3 python3-wheel  python2-dev python-pip python-setuptools

WORKDIR /radar
ADD . /radar
RUN bash build.sh


FROM ubuntu:latest AS runtime-image
RUN apt-get update
RUN apt-get install -y --no-install-recommends python2 python-setuptools

COPY --from=compile-image /usr/local/bin/lfasta /usr/local/bin/lfasta
COPY --from=compile-image /usr/local/bin/radar.py /usr/local/bin/radar.py
COPY --from=compile-image /usr/local/lib/python2.7/dist-packages/Radar-1.3-py2.7-linux-x86_64.egg /usr/local/lib/python2.7/dist-packages/Radar-1.3-py2.7-linux-x86_64.egg

CMD ['/usr/local/bin/radar.py']

ENTRYPOINT ['/usr/bin/python2']