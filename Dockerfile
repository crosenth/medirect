# golob/medirect
#
# VERSION               0.8.0_BCW_0.2.0

FROM      ubuntu:16.04
RUN mkdir /fh
RUN mkdir /app
RUN apt-get update && apt-get install -y \
    python3-dev \
    python3-pip \
    wget \
    perl \
    cpanminus \
    libssl-dev

RUN pip3 install awscli boto3 bucket_command_wrapper==0.2.0

RUN mkdir /src
WORKDIR /src
RUN wget ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/edirect.tar.gz && tar xzvf edirect.tar.gz
RUN cp /src/edirect/* /usr/local/bin/ 
RUN rm -rf /src/edirect/ && rm edirect.tar.gz
RUN cpanm LWP::Protocol::https

RUN mkdir /logs && mkdir /records

ADD medirect/ftract.py /usr/local/bin/ftract
ADD medirect/mefetch.py /usr/local/bin/mefetch
ADD utils/accession_version.py /usr/local/bin/accession_version
ADD utils/ncbi_get_nt_accessions_for_query.py /usr/local/bin/ncbi_get_nt_accessions_for_query

RUN chmod +x /usr/local/bin/*

RUN ln -s /usr/bin/python3 /usr/bin/python
RUN mkdir -p /mnt/inputs/file && mkdir -p /mnt/outputs/file && mkdir /scratch
WORKDIR /root/
