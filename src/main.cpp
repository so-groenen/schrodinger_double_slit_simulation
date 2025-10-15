#include <iostream>
#include <cmath>
#include <print>
#include "raylib.h"
#include "schrodinger_equation.hpp"
#include "interferometer.hpp"

constexpr uint32_t WINDOW_HEIGHT      = 576;
constexpr uint32_t WINDOW_WIDTH       = 1024;
constexpr uint32_t TITLE_FONT_SIZE    = 24;
constexpr uint32_t SUBTITLE_FONT_SIZE = 16;
constexpr uint8_t MAX_COLOR           = 255;

 
enum class SCENE
{
    TITLESCREEN,
    SIMULATION,
};
SCENE current_scene = SCENE::TITLESCREEN;
void show_title_screen(const char* title, const char* subtitle);


int main()
{
    // ===============================================//
    //                                                //
    //           Define system parameters             //    
    //                                                //
    // ===============================================//

    constexpr Vector2 L  = {.x = 6.f,     .y = 6.f};     // system size
    constexpr Vector2 dr = {.x = 0.05f,   .y = 0.05f};   // step size
    constexpr Vector2 r0 = {.x = L.x/5.f, .y = L.y/2.f}; // initial position of the wf

    SchodingerEquation schrodinger{L, dr, r0};
    size_t Nx = schrodinger.Nx;
    size_t Ny = schrodinger.Ny;

    Interferometer double_slit{Nx, Ny};

  
    // ================================================================== //
    //                                                                    //
    //     Set the parameter (in pixel) for the Young interferometer      //    
    //                                                                    //
    // ================================================================== //

    constexpr size_t thickness = 2;
    constexpr size_t width     = 21;
    constexpr size_t opening   = 4;
    size_t height    = (Ny)/2 - width/2 - opening;
    
    double_slit.set_param(thickness, opening, width, height);
 
   
    // ===============================================//
    //                                                //
    //            Raylib & Draw parameters            //    
    //                                                //
    // ===============================================//
    InitWindow(WINDOW_WIDTH, WINDOW_HEIGHT, "Double Slit Experiment");

    float system_width = GetScreenWidth()*0.9f;
    float system_heigh = GetScreenHeight()*0.9f;

    // Size in pixel of the tiles on which we draw the wavefunction
    uint32_t dxR  = ( system_width/(float)Nx );
    uint32_t dyR  = ( system_heigh/(float)Ny ) ;
 
    // parameters for drawing
    uint32_t xStart = (GetScreenWidth()  - Nx*dxR )/2;
    uint32_t yStart = (GetScreenHeight() - Ny*dyR )/2;
    uint32_t midX   = Nx/2;
    uint32_t midY   = Ny/2;
    uint32_t top    = 0;
    uint32_t bottom = Ny-1;

    // Quantum box
    Rectangle quantum_box{(float)xStart, (float)yStart, (Nx)*(float)dxR , (Ny)*(float)dyR};

    // Create small tile to draw on; load it in the GPU for faster drawing.
    RenderTexture2D tile = LoadRenderTexture(dxR, dyR);
    BeginTextureMode(tile);
        ClearBackground(WHITE);
        DrawRectangle(0, 0, dxR, dyR, WHITE);
    EndTextureMode();
    
    // Simulation variables:
    float value{};                                         // local value of the wavefunction .
    float vmax{};                                          // max value of the wavefunction at time t-1
    float vmax_new{};                                      // max value of the wavefunction at time t
    float alpha{};                                         // value of transparency for the tiles
    Color tile_color{MAX_COLOR, 100, 100, (uint8_t)alpha};  // tile color, initiliazed with alpha = 0
    uint32_t iy{}, jx{};


    // ===============================================//
    //                                                //
    //              START THE SIMULATION              //    
    //                                                //
    // ===============================================//

    schrodinger.interact(double_slit);
    schrodinger.evolve();
    vmax = schrodinger.get_max_amplitude();

    size_t Ntotal = (Ny-2)*(Nx-2); // total number of non zero elements

    while (!WindowShouldClose())
    {
        if(IsKeyPressed(KEY_ENTER) || IsKeyPressed(KEY_SPACE))
        {
            current_scene = SCENE::SIMULATION;
        }
        if(current_scene == SCENE::TITLESCREEN)
        {
            static const char* title    = "Young Experiment: Solving the 2D Schrödinger equation";
            static const char* subtitle = "Press SPACE or ENTER";
            show_title_screen(title, subtitle);
        }    
        if(current_scene == SCENE::SIMULATION)
        {
            schrodinger.interact(double_slit);


            BeginDrawing();
                ClearBackground(Color{30, 30, 30, MAX_COLOR});
                for (size_t k = 0; k < Ntotal; k++)
                {
                    value       = schrodinger.get_wf_modulus(k);
                    alpha       = (MAX_COLOR*value/vmax);
                    tile_color.a = (alpha > MAX_COLOR)? MAX_COLOR : (uint8_t)alpha;

                    if (value > vmax_new) vmax_new = value;

                    iy = k%(Ny-2);          
                    jx = k/(Ny-2);

                    uint32_t slit_indx_x = midX - double_slit.thickness/2  + jx;
                    if( (iy < double_slit.width ) && (jx < double_slit.thickness) )
                    {
                        uint32_t slit_indx_y = midY - double_slit.width/2  + iy;
                        DrawTexture(tile.texture, xStart + slit_indx_x*dxR, yStart + slit_indx_y*dyR, RAYWHITE);
                    }
                    if( (iy < double_slit.height ) && (jx < double_slit.thickness))
                    {
                        uint32_t top_part_y    = top    + iy;                           
                        uint32_t bottom_part_y = bottom - iy;                           
                        DrawTexture(tile.texture, xStart + slit_indx_x*dxR, yStart +  top_part_y   *dyR, RAYWHITE);
                        DrawTexture(tile.texture, xStart + slit_indx_x*dxR, yStart +  bottom_part_y*dyR, RAYWHITE);
                    }
                    DrawTexture(tile.texture, xStart + (jx+1)*dxR, yStart + (iy+1)*dyR, tile_color);
                }               
                DrawRectangleLinesEx(quantum_box, 4, WHITE);
                DrawFPS(xStart, yStart*0.4);
            EndDrawing();

            vmax     = vmax_new;
            vmax_new = 0;
            
            schrodinger.evolve();
        }
    }

    UnloadRenderTexture(tile);
    CloseWindow();    
    return EXIT_SUCCESS;
}

 
void show_title_screen(const char* title, const char* subtitle)
{
    int title_size       = MeasureText(title, TITLE_FONT_SIZE);
    int subtitle_size    = MeasureText(subtitle, SUBTITLE_FONT_SIZE);

    BeginDrawing();
        ClearBackground(BLACK);
        DrawText(title,    (GetScreenWidth()-title_size   )/2, GetScreenHeight()*0.4,    TITLE_FONT_SIZE, WHITE);
        DrawText(subtitle, (GetScreenWidth()-subtitle_size)/2, GetScreenHeight()*0.5, SUBTITLE_FONT_SIZE, WHITE);
    EndDrawing(); 
}
