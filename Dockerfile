# Use an official Python slim image. 
# python:3.10-slim satisfies your ">=3.8" requirement and is lightweight.
FROM python:3.10-slim

# Set environment variables for best practices in containers
# 1. Prevents Python from writing .pyc files
ENV PYTHONDONTWRITEBYTECODE=1
# 2. Ensures Python output is sent straight to the terminal (essential for logs)
ENV PYTHONUNBUFFERED=1

# Set the working directory inside the container
WORKDIR /app

RUN apt-get update && \
    apt-get install -y --no-install-recommends procps && \
    rm -rf /var/lib/apt/lists/*
    
# Upgrade pip
RUN pip install --no-cache-dir --upgrade pip

# Copy the file that defines your project and dependencies
COPY pyproject.toml .

# Copy your application's source code, as defined in your toml
# [tool.hatchling.build.targets.wheel] packages = ["src/ReAnnota"]
COPY src ./src

# Copy the README, as it's referenced in your toml
COPY README.md .

# Install the project.
# This one command will:
# 1. Read pyproject.toml
# 2. Install the build backend (hatchling)
# 3. Install all runtime dependencies from [project.dependencies]
# 4. Build and install your "ReAnnota" package
# It will *not* install your [project.optional-dependencies]dev, which is perfect.
RUN pip install --no-cache-dir .

# --- Security Best Practice ---
# Create a non-root user to run the application
RUN useradd -m -u 1000 appuser
USER appuser

# --- Run the App ---
# Set the entrypoint to "reannot"
# This is the script name defined in your [project.scripts]
# This makes your container behave like the "reannot" executable.
# You can run `docker run my-image --version` or `docker run my-image --help`
#ENTRYPOINT ["reannota"]

# The default command is empty, so the entrypoint runs.
# Users can pass arguments to `docker run` which will be passed to "reannot"
#CMD []