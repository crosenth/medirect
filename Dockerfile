FROM      alpine:3.7
RUN apk add --no-cache  bash \
                        python3 \
                        python3-dev \
                        perl \
                        build-base \
                        openssl \
                        perl-dev \
                        perl-crypt-ssleay==0.72-r7 \
                        perl-mozilla-ca==20160104-r0 \
                        perl-lwp-protocol-https==6.06-r1 \
                        perl-test-requiresinternet==0.05-r0 \
                        expat-dev

RUN ln -s /usr/bin/python3 /usr/local/bin/python
RUN pip3 install pip --upgrade && pip install wheel
RUN pip3 install wheel \
        awscli>=1.15.14 \
        boto3>=1.7.14 \
        numpy>=1.14.2 \
        bucket_command_wrapper==0.3.1 \
        biopython>=1.68

WORKDIR /src
RUN wget ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/edirect.tar.gz && tar xzvf edirect.tar.gz
RUN cp /src/edirect/* /usr/local/bin/ 
RUN rm -rf /src/edirect/ && rm edirect.tar.gz

RUN wget -O - http://cpanmin.us | perl - App::cpanminus
RUN cpanm --installdeps XML::Simple && cpanm XML::Simple

RUN mkdir -p /src/medirect
WORKDIR /src/
COPY setup.py /src/setup.py
COPY medirect/__init__.py /src/medirect/__init__.py
COPY medirect/mefetch.py /src/medirect/mefetch.py
COPY medirect/ftract.py /src/medirect/ftract.py
RUN python3 setup.py install

ADD utils/accession_version.py /usr/local/bin/accession_version
ADD utils/ncbi_get_nt_accessions_for_query.py /usr/local/bin/ncbi_get_nt_accessions_for_query
ADD utils/extract_genbank.py /usr/local/bin/extract_genbank
RUN chmod +x /usr/local/bin/*

RUN ln -s /usr/bin/python3 /usr/bin/python

RUN mkdir -p /fh && mkdir -p /app && mkdir -p /src
RUN mkdir -p /mnt/inputs/file && mkdir -p /mnt/outputs/file && mkdir /scratch
RUN mkdir /logs && mkdir /records

WORKDIR /working