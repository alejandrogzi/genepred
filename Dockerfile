# ---------- Build Stage ----------
FROM rust:1.93.0-bookworm AS builder

WORKDIR /app

COPY genepred/Cargo.toml genepred/Cargo.lock ./
COPY genepred/src ./src

RUN cargo build --release --all-features --bin genepred --locked && \
    strip target/release/genepred

# ---------- Runtime Stage ----------
FROM debian:bookworm-slim

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    ca-certificates \
    procps \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /app/target/release/genepred /usr/local/bin/genepred

# Set up non-root user
RUN useradd -m -u 1000 cuser && \
    chmod +x /usr/local/bin/genepred

USER cuser
WORKDIR /data

RUN genepred --version
RUN genepred lint --help

CMD ["bash"]
