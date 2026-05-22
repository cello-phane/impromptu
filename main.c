#include <SDL.h>

#include "engine.h"
#include "obj_parse.h"

int main(int argc, char *argv[]) {
    //struct Engine *e = Engine_create(3840 / 2, 2160 / 2);

    struct Engine *e = Engine_create(1920, 1080);
    Engine_run(e);
    Engine_destroy(e);

    return 0;
}
