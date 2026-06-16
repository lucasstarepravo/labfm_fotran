FROM torchfort:gpu_cpu_only

RUN apt-get update && apt-get install -y libopenblas-dev

WORKDIR /app