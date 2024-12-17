FROM bojobo/sas:21.0.0 AS base

USER 0

RUN apt-get update && apt-get upgrade -y && apt-get dist-upgrade -y \
    && apt-get install -y --no-install-recommends git \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages
RUN uv pip install beautifultable h5py loguru lxml pydantic pypdf python-dotenv yt pre-commit \
    && uv cache clean

FROM base

# Add SIXTE
COPY --from=bojobo/sixte:3.0 --chown=heasoft:heasoft /opt/simput /opt/simput
ENV SIMPUT=/opt/simput \
    SIXTE=/opt/simput \
    PATH=/opt/simput/bin:${PATH} \
    PFILES=${PFILES}:/opt/simput/share/sixte/pfiles:/opt/simput/share/simput/pfiles \
    LD_LIBRARY_PATH=/opt/simput/lib:${LD_LIBRARY_PATH}

USER heasoft
