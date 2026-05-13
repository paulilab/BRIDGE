#!/bin/sh
# Fix write permissions on the mounted database so the shiny user can write
# the storr cache tables into it. Runs as root before dropping privileges.
if [ -f /srv/data/database.db ]; then
    chmod 666 /srv/data/database.db
fi

exec su shiny -s /bin/sh -c \
    "R -e \"options(shiny.host='0.0.0.0', shiny.port=3838); shiny::runApp('/srv/shiny-server/app.R')\""