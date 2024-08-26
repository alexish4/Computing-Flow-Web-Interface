# wsgi.py

from application import app  # Import the Flask app instance from application.py

if __name__ == "__main__":
    app.run()
# import sys
# import os

# # Add the directory containing 'application.py' to the Python path
# sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# from application import app  # Import the Flask app instance
