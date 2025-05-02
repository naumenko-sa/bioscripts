#!/bin/sh

# download app results - dragen output only

# list appsessions for a project
bs list appsessions --project-name=[project_name]

# list datasets in appsession
bs list datasets --filter-field AppSession.Id --filter-term [session_id]

bs download dataset -i [ds.id] -o .

