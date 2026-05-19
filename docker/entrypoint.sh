#!/bin/sh
set -eu

APP_CMD="options(shiny.host='0.0.0.0', shiny.port=3838); shiny::runApp('/srv/shiny-server/app.R')"

# If we are root, keep the previous behavior:
# 1) fix bind-mounted DB permissions for write access
# 2) drop privileges to shiny
if [ "$(id -u)" -eq 0 ]; then
    if [ -f /srv/data/database.db ]; then
        chmod 666 /srv/data/database.db || true
    fi

    if id shiny >/dev/null 2>&1; then
        exec su shiny -s /bin/sh -c "R -e \"${APP_CMD}\""
    fi
fi

# OpenShift typically runs with an arbitrary non-root UID and no password entry,
# so `su` is not usable there. Run directly as current user instead.
exec R -e "${APP_CMD}"