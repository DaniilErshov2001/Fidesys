#ifndef CONSTRAINT_H
#define CONSTRAINT_H

struct Constraint
{
    enum Type
    {
        NONE = 0,
        UX = 1,
        UY = 2,
        UXUY = 3,
    };
    int node;
    Type type;
};

#endif // CONSTRAINT_H
