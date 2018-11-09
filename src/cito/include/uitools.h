/*  Copyright © 2018, Roboti LLC

    This file is licensed under the MuJoCo Resource License (the "License").
    You may not use this file except in compliance with the License.
    You may obtain a copy of the License at

        https://www.roboti.us/resourcelicense.txt
*/


#pragma once


#include "mujoco.h"
#include "glfw3.h"


// this is a C-API
#if defined(__cplusplus)
extern "C"
{
#endif


// User-supplied callback function types.
typedef void (*uiEventFn)(mjuiState* state);
typedef void (*uiLayoutFn)(mjuiState* state);

// Container for GLFW window pointer.
struct _uiUserPointer
{
    mjuiState* state;
    uiEventFn uiEvent;
    uiLayoutFn uiLayout;
    double buffer2window;
};
typedef struct _uiUserPointer uiUserPointer;

// Set internal and user-supplied UI callbacks in GLFW window.
void uiSetCallback(GLFWwindow* wnd, mjuiState* state, 
                   uiEventFn uiEvent, uiLayoutFn uiLayout);

// Clear UI callbacks in GLFW window.
void uiClearCallback(GLFWwindow* wnd);

// Compute suitable font scale.
int uiFontScale(GLFWwindow* wnd);

// Modify UI structure.
void uiModify(GLFWwindow* wnd, mjUI* ui, mjuiState* state, mjrContext* con);


#if defined(__cplusplus)
}
#endif
