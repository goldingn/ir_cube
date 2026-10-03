#!/bin/bash
# Start sshd if PUBLIC_KEY is set (RunPod sets it from the account's SSH keys),
# then run the command (by default, sleep).
set -e

if [ -n "${PUBLIC_KEY}" ]; then
  mkdir -p /root/.ssh
  chmod 700 /root/.ssh
  printf '%s\n' "${PUBLIC_KEY}" > /root/.ssh/authorized_keys
  chmod 600 /root/.ssh/authorized_keys

  # ssh sessions do not inherit the container's environment: pass it on through
  # /etc/environment (read by PAM for every session), so that RETICULATE_PYTHON,
  # the thread counts and anything set in the pod template apply there too
  printenv \
    | grep -vE '^(HOME|HOSTNAME|PWD|OLDPWD|SHLVL|TERM|PUBLIC_KEY|_)=' \
    | grep -E '^[A-Za-z_][A-Za-z0-9_]*=[^"]*$' \
    | sed -E 's/^([^=]+)=(.*)$/\1="\2"/' \
    > /etc/environment

  ssh-keygen -A
  mkdir -p /run/sshd
  /usr/sbin/sshd
  echo "sshd started"
else
  echo "PUBLIC_KEY is not set: sshd not started"
fi

exec "$@"
