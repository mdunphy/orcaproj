#ifndef ORCA_TRIPOLAR_INIT_H
#define ORCA_TRIPOLAR_INIT_H

// Initializes the library (idempotent)
// May be called explicitly, but API functions will initialize automatically if needed
int orca_init(void);

// Optional: set data directory
// Must be called before the first call to orca_init() or any API function
// Calling this after initialization has no effect
void orca_set_data_dir(const char *dir);

#endif /* ORCA_TRIPOLAR_INIT_H */
