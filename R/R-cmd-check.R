# fix the note:
# no visible binding for global variable '.'
utils::globalVariables(".")


# fix the note:
# no visible global function definition for 'between, site, start, end'
utils::globalVariables(c("between", "site", "start", "end"))
