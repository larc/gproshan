#ifndef GLFW_KEYS_H
#define GLFW_KEYS_H


#include <gproshan/viewer/include_opengl.h>

#include <unordered_map>


// geometry processing and shape analysis framework
namespace gproshan {


const std::unordered_map<int, const char *> glfw_key_name = {
															{GLFW_KEY_UNKNOWN, ""},//-1
															{GLFW_KEY_SPACE, "SPACE"},//32
															{GLFW_KEY_APOSTROPHE, "'"},//39 /* ' */
															{GLFW_KEY_COMMA, ","},//44 /* , */
															{GLFW_KEY_MINUS, "-"},//45 /* - */
															{GLFW_KEY_PERIOD, "."},//46 /* . */
															{GLFW_KEY_SLASH, "/"},//47 /* / */
															{GLFW_KEY_0, "0"},//48
															{GLFW_KEY_1, "1"},//49
															{GLFW_KEY_2, "2"},//50
															{GLFW_KEY_3, "3"},//51
															{GLFW_KEY_4, "4"},//52
															{GLFW_KEY_5, "5"},//53
															{GLFW_KEY_6, "6"},//54
															{GLFW_KEY_7, "7"},//55
															{GLFW_KEY_8, "8"},//56
															{GLFW_KEY_9, "9"},//57
															{GLFW_KEY_SEMICOLON, ";"},//59 /* ; */
															{GLFW_KEY_EQUAL, "="},//61 /* = */
															{GLFW_KEY_A, "A"},//65
															{GLFW_KEY_B, "B"},//66
															{GLFW_KEY_C, "C"},//67
															{GLFW_KEY_D, "D"},//68
															{GLFW_KEY_E, "E"},//69
															{GLFW_KEY_F, "F"},//70
															{GLFW_KEY_G, "G"},//71
															{GLFW_KEY_H, "H"},//72
															{GLFW_KEY_I, "I"},//73
															{GLFW_KEY_J, "J"},//74
															{GLFW_KEY_K, "K"},//75
															{GLFW_KEY_L, "L"},//76
															{GLFW_KEY_M, "M"},//77
															{GLFW_KEY_N, "N"},//78
															{GLFW_KEY_O, "O"},//79
															{GLFW_KEY_P, "P"},//80
															{GLFW_KEY_Q, "Q"},//81
															{GLFW_KEY_R, "R"},//82
															{GLFW_KEY_S, "S"},//83
															{GLFW_KEY_T, "T"},//84
															{GLFW_KEY_U, "U"},//85
															{GLFW_KEY_V, "V"},//86
															{GLFW_KEY_W, "W"},//87
															{GLFW_KEY_X, "X"},//88
															{GLFW_KEY_Y, "Y"},//89
															{GLFW_KEY_Z, "Z"},//90
															{GLFW_KEY_LEFT_BRACKET, "["},//91 /* [ */
															{GLFW_KEY_BACKSLASH, "\\"},//92 /* \ */
															{GLFW_KEY_RIGHT_BRACKET, "]"},//93 /* ] */
															{GLFW_KEY_GRAVE_ACCENT, "`"},//96 /* ` */
															{GLFW_KEY_WORLD_1, ""},//161 /* non-US #1 */
															{GLFW_KEY_WORLD_2, ""},//162 /* non-US #2 */
															{GLFW_KEY_ESCAPE, "ESCAPE"},//256
															{GLFW_KEY_ENTER, "ENTER"},//257
															{GLFW_KEY_TAB, "TAB"},//258
															{GLFW_KEY_BACKSPACE, "BACKSPACE"},//259
															{GLFW_KEY_INSERT, "INSERT"},//260
															{GLFW_KEY_DELETE, "DELETE"},//261
															{GLFW_KEY_RIGHT, "RIGHT"},//262
															{GLFW_KEY_LEFT, "LEFT"},//263
															{GLFW_KEY_DOWN, "DOWN"},//264
															{GLFW_KEY_UP, "UP"},//265
															{GLFW_KEY_PAGE_UP, ""},//266
															{GLFW_KEY_PAGE_DOWN, ""},//267
															{GLFW_KEY_HOME, ""},//268
															{GLFW_KEY_END, ""},//269
															{GLFW_KEY_CAPS_LOCK, ""},//280
															{GLFW_KEY_SCROLL_LOCK, ""},//281
															{GLFW_KEY_NUM_LOCK, ""},//282
															{GLFW_KEY_PRINT_SCREEN, ""},//283
															{GLFW_KEY_PAUSE, ""},//284
															{GLFW_KEY_F1, "F1"},//290
															{GLFW_KEY_F2, "F2"},//291
															{GLFW_KEY_F3, "F3"},//292
															{GLFW_KEY_F4, "F4"},//293
															{GLFW_KEY_F5, "F5"},//294
															{GLFW_KEY_F6, "F6"},//295
															{GLFW_KEY_F7, "F7"},//296
															{GLFW_KEY_F8, "F8"},//297
															{GLFW_KEY_F9, "F9"},//298
															{GLFW_KEY_F10, "F10"},//299
															{GLFW_KEY_F11, "F11"},//300
															{GLFW_KEY_F12, "F12"},//301
															{GLFW_KEY_F13, ""},//302
															{GLFW_KEY_F14, ""},//303
															{GLFW_KEY_F15, ""},//304
															{GLFW_KEY_F16, ""},//305
															{GLFW_KEY_F17, ""},//306
															{GLFW_KEY_F18, ""},//307
															{GLFW_KEY_F19, ""},//308
															{GLFW_KEY_F20, ""},//309
															{GLFW_KEY_F21, ""},//310
															{GLFW_KEY_F22, ""},//311
															{GLFW_KEY_F23, ""},//312
															{GLFW_KEY_F24, ""},//313
															{GLFW_KEY_F25, ""},//314
															{GLFW_KEY_KP_0, ""},//320
															{GLFW_KEY_KP_1, ""},//321
															{GLFW_KEY_KP_2, ""},//322
															{GLFW_KEY_KP_3, ""},//323
															{GLFW_KEY_KP_4, ""},//324
															{GLFW_KEY_KP_5, ""},//325
															{GLFW_KEY_KP_6, ""},//326
															{GLFW_KEY_KP_7, ""},//327
															{GLFW_KEY_KP_8, ""},//328
															{GLFW_KEY_KP_9, ""},//329
															{GLFW_KEY_KP_DECIMAL, ""},//330
															{GLFW_KEY_KP_DIVIDE, ""},//331
															{GLFW_KEY_KP_MULTIPLY, ""},//332
															{GLFW_KEY_KP_SUBTRACT, ""},//333
															{GLFW_KEY_KP_ADD, ""},//334
															{GLFW_KEY_KP_ENTER, ""},//335
															{GLFW_KEY_KP_EQUAL, ""},//336
															{GLFW_KEY_LEFT_SHIFT, ""},//340
															{GLFW_KEY_LEFT_CONTROL, ""},//341
															{GLFW_KEY_LEFT_ALT, ""},//342
															{GLFW_KEY_LEFT_SUPER, ""},//343
															{GLFW_KEY_RIGHT_SHIFT, ""},//344
															{GLFW_KEY_RIGHT_CONTROL, ""},//345
															{GLFW_KEY_RIGHT_ALT, ""},//346
															{GLFW_KEY_RIGHT_SUPER, ""},//347
															{GLFW_KEY_MENU, ""},//348
															{GLFW_KEY_LAST, ""},//GLFW_KEY_MENU
															};


} // namespace gproshan

#endif // GLFW_KEYS_H

