#!/bin/sh

set -e

download() {
  if command -v curl > /dev/null 2>&1; then
    curl -fsSL "$1"
  else
    wget -qO- "$1"
  fi
}

APP_NAME=nf-test
GITHUB_ORG=askimed
GITHUB_REPO=nf-test

if [ -n "$1" ]; then
  VERSION="$1"
else
  GITHUB_LATEST_RELEASE_URL=https://api.github.com/repos/${GITHUB_ORG}/${GITHUB_REPO}/releases/latest
  VERSION_JSON="$(download ${GITHUB_LATEST_RELEASE_URL})"
  VERSION="$(printf '%s' "${VERSION_JSON}" |  awk -F '"' '/tag_name/{print $4}')"
  #remove v prefix
  VERSION="${VERSION:1}"
fi

GITHUB_REPO_URL=https://github.com/${GITHUB_ORG}/${GITHUB_REPO}
GITHUB_RELEASE_URL=${GITHUB_REPO_URL}/releases/download/v${VERSION}/${APP_NAME}-${VERSION}.tar.gz

# download and extract tar.gz file
echo "Downloading ${APP_NAME} ${VERSION} from ${GITHUB_RELEASE_URL}..."
download ${GITHUB_RELEASE_URL} | tar -xz

# move jar file to .nf-test folder
APP_HOME=${HOME}/.${APP_NAME}
mkdir -p ${APP_HOME}
mv -f ${APP_NAME}.jar ${APP_HOME}/${APP_NAME}.jar

echo ""
echo "Done."
