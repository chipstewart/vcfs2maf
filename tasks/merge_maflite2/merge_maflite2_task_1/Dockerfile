
FROM ubuntu:16.04

RUN apt-get update && apt-get install -y python sudo dstat wget zlib1g zlib1g-dev

RUN apt-get update \
	&& apt-get install -y --no-install-recommends ca-certificates \
	&& rm -rf /var/lib/apt/lists/*
	&& rm -rf /tmp/*

ENV JULIA_PATH /usr/local/julia
#ENV JULIA_VERSION 0.5.0
ENV JULIA_VERSION 1.0.5

RUN mkdir $JULIA_PATH \
	&& apt-get update && apt-get install -y curl \
	&& curl -sSL "https://julialang.s3.amazonaws.com/bin/linux/x64/${JULIA_VERSION%[.-]*}/julia-${JULIA_VERSION}-linux-x86_64.tar.gz" -o julia.tar.gz \
	&& curl -sSL "https://julialang.s3.amazonaws.com/bin/linux/x64/${JULIA_VERSION%[.-]*}/julia-${JULIA_VERSION}-linux-x86_64.tar.gz.asc" -o julia.tar.gz.asc \
	&& export GNUPGHOME="$(mktemp -d)" \
# http://julialang.org/juliareleases.asc
# Julia (Binary signing key) <buildbot@julialang.org>
	&& gpg --keyserver ha.pool.sks-keyservers.net --recv-keys 3673DF529D9049477F76B37566E3C7DC03D6E495 \
	&& gpg --batch --verify julia.tar.gz.asc julia.tar.gz \
	&& rm -r "$GNUPGHOME" julia.tar.gz.asc \
	&& tar -xzf julia.tar.gz -C $JULIA_PATH --strip-components 1 \
	&& rm -rf /var/lib/apt/lists/* julia.tar.gz* \
	&& rm -rf /tmp/*


ENV PATH $JULIA_PATH/bin:$PATH

RUN julia -e 'Pkg.add("DataFrames")'
RUN julia -e 'Pkg.build("DataFrames")'
RUN julia -e 'Pkg.add("Missings")'
RUN julia -e 'Pkg.build("Missings")'
RUN julia -e 'Pkg.add("Printf")'
RUN julia -e 'Pkg.build("Printf")'
RUN julia -e 'Pkg.add("DelimitedFiles")'
RUN julia -e 'Pkg.build("DelimitedFiles")'
RUN julia -e 'Pkg.update()'

#copy contents of tasks/<taskname>/build/src on the build host into /opt/src on the docker
COPY src/ /src

WORKDIR /src

