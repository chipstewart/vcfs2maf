FROM julia:1.0.5

#RUN julia -e 'using Pkg; Pkg.REPLMode.pkgstr("add CSV        ;precompile");using CSV'
#RUN julia -e 'using Pkg; Pkg.REPLMode.pkgstr("add DataFrames ;precompile");using DataFrames'
#RUN julia -e 'using Pkg; Pkg.REPLMode.pkgstr("add Missings   ;precompile");using Missings'
#RUN julia -e 'using Pkg; Pkg.REPLMode.pkgstr("add Printf     ;precompile");using Printf'
#RUN julia -e 'using Pkg; Pkg.REPLMode.pkgstr("add DelimitedFiles     ;precompile");using DelimitedFiles'
#RUN julia -e 'Pkg.update()'

RUN julia -e 'import Pkg; Pkg.add("DataFrames"); using DataFrames'
RUN julia -e 'import Pkg; Pkg.add("CSV"); using CSV'
RUN julia -e 'import Pkg; Pkg.add("Missings"); using Missings'
RUN julia -e 'import Pkg; Pkg.add("Printf"); using Printf'
RUN julia -e 'import Pkg; Pkg.add("DelimitedFiles"); using DelimitedFiles'

RUN apt-get update \
	&& apt-get install -y --no-install-recommends ca-certificates \
	&& rm -rf /var/lib/apt/lists/* \
	&& rm -rf /tmp/*

RUN apt-get update && apt-get install -y sudo wget zlib1g zlib1g-dev

COPY src/ /src
