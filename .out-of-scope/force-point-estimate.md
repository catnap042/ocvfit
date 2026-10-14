# Reporting a single number for a non-identifiable mode

Requests for an option that always prints a point value, for example "just give me
LAM_PE for my LFP cell", are out of scope.

## Why this is out of scope

The value in parentheses is the best fit, but when the data do not determine the
mode, other values fit just as well. Printing it without the interval would present
noise as a result. The JSON report already contains the raw best-fit value for anyone
who needs it, next to the verdict.

To make a mode identifiable, change the experiment rather than the report: discharge
deeper so the steep end of the flat electrode is reached, or add a harvested-electrode
half-cell.
