echo "Starting backup process."
# Create a backup directory if it doesn't exist
BACKUP_DIR="backup_$(date +%Y%m%d_%H%M%S)"
mkdir "$BACKUP_DIR"

# Backup include/externs
if [ -d "include/externs" ]; then
  cp -r include/externs "$BACKUP_DIR/"
  echo "Backed up include/externs to $BACKUP_DIR/externs"
else
  echo "Warning: include/externs directory not found."
fi

# Backup src
if [ -d "src" ]; then
  cp -r src "$BACKUP_DIR/"
  echo "Backed up src to $BACKUP_DIR/src"
else
  echo "Warning: src directory not found."
fi

echo "Backup process complete."
