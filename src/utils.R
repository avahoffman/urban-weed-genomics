# Utility functions for the project
#

right_string <-
  function(x, n) {
    # Pulls n number of characters from the end of a string
    # 
    # Args: x: string in question (e.g., file path)
    # n: number of characters back to go
    # 
    substring(x, nchar(x) - n + 1)
  }
