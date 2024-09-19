################################################################################
#                                                                              #
#  Function for custom styling of tables in R Mardown                          #
#  Apply consistent styling for HTML output                                    #
################################################################################

style_kable <- function(data, caption = "") {
  # Format the caption with a custom class for full-width styling
  # styled_caption <- paste0('<span class="kable-caption">', caption, '</span>')
  styled_caption <- paste0('<div style="font-size: 14pt; color: black;
                           font-weight: bold; text-align: center;">',
                           caption, '</div>')
  # Create the table with the styled caption
  kable(data, caption = styled_caption, escape = FALSE) %>%
    # Apply styling to the table
    kable_styling(full_width = FALSE, position = "left", font_size = 12, 
                  bootstrap_options = c("striped", "hover")) %>%
    # Make the header row bold and set background to white
    row_spec(0, bold = TRUE, background = "white") %>%
    # Set background to white for the ID column, make it bold, and set width to 9em
    column_spec(1, width = "9em", bold = TRUE, background = "white",
                extra_css = "white-space: nowrap;") %>%
    # Set background to white for other columns
    column_spec(2:ncol(data), background = "white")
}

# Example usage (you would replace this with your actual data and caption)
# style_kable(head(Data), caption = "Summary of the first few rows of the data")

################################################################################
#                                                                              #
#  Function to record messages in a format for R Markdown                      #
#                                                                              #
################################################################################

logMessage <- function(log_messages, echo = TRUE) {
  
  # Convert messages to a character vector if it's not already
  if (!is.vector(log_messages) || is.numeric(log_messages)) {
    log_messages <- as.character(log_messages)
  }
  
  # Echo to console/HTML
  if (echo) {
    if (knitr::is_html_output()) {
      # Use knitr::asis_output to ensure HTML is not escaped
      knitr::asis_output(
        paste0(
          "<div style='border: 1px solid grey; background-color: white; color: black; padding: 10px; ",
          "margin: 10px 0; border-radius: 10px;'>",  # Added border-radius for rounded corners
          gsub("\n", "<br>", log_messages),
          "</div>"
        )
      )
    } else {
      # Use regular line breaks otherwise
      cat(log_messages, sep = "\n")
    }
  }
}

################################################################################
#                                                                              #
# Find the Google Chrome or Chromium executable based on the operating system  #
# The webshot2 package, needed to print the Sankey plot, requires Chrome       #
################################################################################

ensureChrome <- function() {
  
  # Install and load the 'chromote' package if not already available
  if (!requireNamespace("chromote", quietly = TRUE)) {
    install.packages("chromote")
  }
  library(chromote)
  
  # Function to find Chrome path on macOS
  find_chrome_mac <- function() {
    chrome_paths <- c(
      "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",  # Standard installation path
      "~/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"  # Alternate installation path
    )
    
    for (path in chrome_paths) {
      if (file.exists(path)) {
        return(normalizePath(path))
      }
    }
    return(NULL)
  }
  
  # Main logic for ensuring Chrome is available
  if (Sys.info()["sysname"] == "Darwin") {
    # macOS specific handling
    chrome_path <- find_chrome_mac()
    if (!is.null(chrome_path)) {
      Sys.setenv(CHROMOTE_CHROME = chrome_path)
      b <- ChromoteSession$new()
      b$close()
      return(TRUE)
    } else {
      message("Google Chrome was not found in the usual paths on macOS. Please install it or specify the path.")
      return(FALSE)
    }
  } else {
    # For other operating systems, rely on Chromote's default behavior
    tryCatch({
      b <- ChromoteSession$new()
      b$close()
      return(TRUE)
    }, error = function(e) {
      message("Google Chrome was not found. Please install it or set the CHROMOTE_CHROME environment variable.")
      return(FALSE)
    })
  }
}
