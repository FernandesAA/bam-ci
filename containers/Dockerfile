# Imagem equivalente ao containers/bam-ci.def, para hosts com Docker.
#   docker build -f containers/Dockerfile -t bam-ci:1.0.0 .
#   docker run --rm bam-ci:1.0.0 --help
# E possivel converter para Apptainer sem root:
#   apptainer build bam-ci.sif docker-daemon://bam-ci:1.0.0
FROM python:3.12-slim

RUN apt-get update \
 && apt-get install -y --no-install-recommends libeccodes0 \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/bam-ci
COPY pyproject.toml README.md LICENSE ./
COPY src ./src
COPY tests ./tests
RUN python -m pip install --no-cache-dir --upgrade pip \
 && python -m pip install --no-cache-dir ".[grib,netcdf,download]" \
 && python -c "import bamci; print('bam-ci', bamci.__version__)"

ENV LC_ALL=C.UTF-8 LANG=C.UTF-8 PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1

ENTRYPOINT ["bamci"]
CMD ["--help"]
