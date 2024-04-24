FROM --platform=amd64 debian:bookworm-slim
RUN export TZ=Etc/UTC
RUN DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get update && \
apt-get -y install tzdata && \
apt-get install -y \
        wget \
        python3-pip \
        perl \
&& apt-get clean \
&& apt-get purge \
&& rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN mkdir /src/
WORKDIR /src

RUN wget -O - http://cpanmin.us | perl - App::cpanminus && \
cpanm --installdeps XML::Simple && cpanm XML::Simple

RUN wget https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/versions/20.9.20231210/edirect-20.9.20231210.tar.gz && \
tar xzvf edirect-20.9.20231210.tar.gz && \
cp -r /src/edirect/* /usr/local/bin/ && \
rm -rf /src/edirect/ && \
rm edirect-20.9.20231210.tar.gz

ADD ./src/ /src/medirect/src
ADD ./pyproject.toml /src/medirect
RUN cd /src/medirect && \
pip3  install .  --break-system-packages

ADD utils/accession_version.py /usr/local/bin/accession_version
ADD utils/ncbi_get_nt_accessions_for_query.py /usr/local/bin/ncbi_get_nt_accessions_for_query
ADD utils/extract_genbank.py /usr/local/bin/extract_genbank
RUN chmod +x /usr/local/bin/*

RUN mkdir /logs && mkdir /records && mkdir /working
RUN ln -s /usr/bin/python3 /usr/local/bin/python

WORKDIR /working