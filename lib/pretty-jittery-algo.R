library(tidyverse)

prettyJitter <- function(df, col1, col2, col1Iterator) {
  accuracy= 100
  # Using Round(accuracy * log10) to make sure 
  # similar value will be horizontally aligned.
  df <- mutate(df, Log10Value=round(log10(df[[col2]]) * accuracy))
  
  new_df = data.frame()
  
  for (x_idx in 1:length(col1Iterator)) {
    x_name = col1Iterator[x_idx]
    
    # Get rows for one x-axis item
    df_rows = filter(df, df[[col1]] == x_name)
    
    new_df_rows = data.frame()
    
    # Get unique y-axis value
    unique_y_count = count(df_rows, Log10Value)
    
    if (nrow(unique_y_count) == 0) {
      next
    }
    for (y_idx in 1:nrow(unique_y_count)) {
      y_value = as.integer(unique_y_count[y_idx,1])
      y_count = as.integer(unique_y_count[y_idx,2])
      df_rows_with_value = filter(df_rows, Log10Value==y_value)
      
      left_cursor = x_idx
      right_cursor = x_idx
      step_length = 0.05
      even_split_length =0.025
      
      # generate start x-axis position for same y-axis value
      if (y_count %% 2 == 0) {
        left_cursor = x_idx - even_split_length
        right_cursor = x_idx + even_split_length
        x_pos = c(left_cursor, right_cursor)
        y_count = y_count - 2
      } else {
        x_pos = c(x_idx)
        y_count = y_count - 1
      }
      
      # generate all x-axis position for same y-axis value
      for(c in 1:y_count) {
        if(y_count < 1) {
          break
        }
        if (c %% 2) {
          left_cursor = left_cursor - step_length
          x_pos = prepend(x_pos, left_cursor)
        } else {
          right_cursor = right_cursor + step_length
          x_pos = append(x_pos, right_cursor)
        }
      }
      
      # add position to each rows for same y-axis value
      for(c in 1:length(x_pos)) {
        p = x_pos[c]
        row = df_rows_with_value[c,]
        row['X'] <- p
        row['Col2'] = 10 ^ (row$Log10Value / accuracy)
        new_df_rows = rbind(new_df_rows, row)
      }
    }
    
    new_df = rbind(new_df, new_df_rows)
  }
  
  return(new_df)
}
