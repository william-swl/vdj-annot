prefix := '#!/usr/bin/bash'
cellid_pattern := '^ *"id": "[a-z0-9\-]*",$'
ignore_jupyter_cellid_content := "grep -v \\'" + cellid_pattern + "\\'"

init:
    # script
    echo '{{prefix}}' > src/ignore_jupyter_cellid.sh
    echo {{ignore_jupyter_cellid_content}} >> src/ignore_jupyter_cellid.sh
    chmod +x src/ignore_jupyter_cellid.sh
    # gitattributes
    echo "*.ipynb filter=igc" > .gitattributes
    # git config
    git config core.attributesfile .gitattributes
    git config filter.igc.clean src/ignore_jupyter_cellid.sh
    git config filter.igc.smudge cat

