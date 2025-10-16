#include <iostream>
#include <cmath>
#include <print>
#include "raylib.h"
#include "interferometer.hpp"
#include "crank_nicolson_builder.hpp"
#include "gaussian_wavefunction_builder.hpp"
#include "schrodinger_equation_builder.hpp"
#include "schrodinger_equation.hpp"

constexpr uint32_t WINDOW_HEIGHT      = 600;
constexpr uint32_t WINDOW_WIDTH       = 1024;
constexpr uint32_t TITLE_FONT_SIZE    = 24;
constexpr uint32_t SUBTITLE_FONT_SIZE = 16;
constexpr uint8_t MAX_COLOR           = 255;
constexpr uint8_t RGB_RED             = MAX_COLOR;
constexpr uint8_t RGB_GREEN           = 100;
constexpr uint8_t RGB_BLUE            = 100;
 
enum class SCENE
{
    TITLESCREEN,
    SIMULATION,
};
SCENE current_scene = SCENE::TITLESCREEN;
void show_title_screen(std::string_view title, std::string_view subtitle);

int main()
{
    // ===============================================//
    //                                                //
    //           Define system parameters             //    
    //                                                //
    // ===============================================//

    constexpr Vector2 L   {.x = 6.f,     .y = 4.f};     // system size
    constexpr Vector2 dr  {.x = 0.04f,   .y = 0.04f};   // step size
    constexpr Vector2 r0  {.x = L.x/5.f, .y = L.y/2.f}; // initial position of the wf

    auto gaussian_wf_builder   = std::make_unique<GaussianWfBuilder>();   
    auto sparse_matrix_builder = std::make_unique<CrankNicolsonBuilder>(); 

    SchodingerEquationBuilder eq_builder {L, dr, r0, std::move(sparse_matrix_builder), std::move(gaussian_wf_builder)};
    SchodingerEquation        schrodinger{eq_builder.build_equation()};

    const size_t Nx           {schrodinger.Nx}; // number of "pixels" (steps) along the x direction
    const size_t Ny           {schrodinger.Ny}; // number of "pixels" (steps) along the y direction
    const size_t Ntotal       {(Ny-2)*(Nx-2)};  // size of the matrix, with the boundary terms removed.
    Interferometer double_slit{Nx, Ny};

  
    // ================================================================== //
    //                                                                    //
    //     Set the parameter (in pixel) for the Young interferometer      //    
    //                                                                    //
    // ================================================================== //

    const size_t thickness = static_cast<size_t>(0.02f*Nx); // thickness of interfero, along the x-axixs, in pixels
    const size_t width     = static_cast<size_t>(0.18f*Ny); // size of the part between the two slits (along the y-axis)
    const size_t opening   = static_cast<size_t>(0.04f*Ny); // size of the two slits (along the y-axis)
    const size_t height{Ny/2 - width/2 - opening};          // "height" from top/bottom until the slit. 
    
    double_slit.set_param(thickness, width, height);
 
   
    // ===============================================//
    //                                                //
    //            Raylib & Draw parameters            //    
    //                                                //
    // ===============================================//
    InitWindow(WINDOW_WIDTH, WINDOW_HEIGHT, "Double Slit Experiment");

    const float system_width {GetScreenWidth()*0.9f};
    const float system_heigh {GetScreenHeight()*0.9f};

    // Size in pixel of the tiles on which we draw the wavefunction
    const size_t tile_dx  = static_cast<size_t>(system_width/Nx);
    const size_t tile_dy  = static_cast<size_t>(system_heigh/Ny);

    // Create small tile to draw on; load it in the GPU for faster drawing.
    RenderTexture2D tile = LoadRenderTexture(tile_dx, tile_dy);
    BeginTextureMode(tile);
        ClearBackground(WHITE);
        DrawRectangle(0, 0, tile_dx, tile_dy, WHITE);
    EndTextureMode();


    // parameters for drawing
    size_t x_start  = static_cast<size_t>((GetScreenWidth()  - Nx*tile_dx )/2);
    size_t y_start  = static_cast<size_t>((GetScreenHeight() - Ny*tile_dy )/2);
    Point start{x_start,y_start};

    // Quantum box
    Rectangle quantum_box{(float)x_start, (float)y_start, Nx*(float)tile_dx , Ny*(float)tile_dy};


    
    // Simulation variables:
    float value{};                                      // local value of the wavefunction .
    float vmax{};                                       // max value of the wavefunction at time t-1
    float vmax_new{};                                   // max value of the wavefunction at time t
    float alpha{};                                      // value of transparency for the tiles
    Color tile_color{RGB_RED, RGB_BLUE, RGB_GREEN, 0};  // tile color, initiliazed with alpha = 0
    size_t iy{}, jx{};


    // ===============================================//
    //                                                //
    //              START THE SIMULATION              //    
    //                                                //
    // ===============================================//

    schrodinger.interact(double_slit);
    schrodinger.evolve();
    vmax = schrodinger.get_max_amplitude();

    SetTargetFPS(60);
    while (!WindowShouldClose())
    {
        if(IsKeyPressed(KEY_ENTER) || IsKeyPressed(KEY_SPACE))
        {
            current_scene = SCENE::SIMULATION;
        }
        
        if(current_scene == SCENE::TITLESCREEN)
        {
            static std::string_view title    = "Young Experiment: Solving the 2D Schrödinger equation";
            static std::string_view subtitle = "Press SPACE or ENTER";
            show_title_screen(title, subtitle);
        }    
        else if(current_scene == SCENE::SIMULATION)
        {
            schrodinger.interact(double_slit);
            BeginDrawing();
                ClearBackground(Color{30, 30, 30, MAX_COLOR});
                for (size_t k = 0; k < Ntotal; k++)
                {
                    value        = schrodinger.get_wf_modulus(k);
                    alpha        = MAX_COLOR*(value/vmax);
                    tile_color.a = (alpha > MAX_COLOR)? MAX_COLOR : static_cast<uint8_t>(alpha);

                    if (value > vmax_new) 
                        vmax_new = value;

                    iy = k%(Ny-2);          
                    jx = k/(Ny-2);
                    double_slit.draw(tile, start, iy, jx);
                    DrawTexture(tile.texture, start.x + (jx+1)*tile_dx, start.y + (iy+1)*tile_dy, tile_color);
                }               
                DrawRectangleLinesEx(quantum_box, 4, WHITE);
                DrawFPS(x_start, y_start*0.4);
            EndDrawing();

            vmax     = vmax_new;
            vmax_new = 0;
            
            schrodinger.evolve();

            if(IsKeyPressed(KEY_BACKSPACE))
            {
                schrodinger.reset();
            }
        }
    }

    UnloadRenderTexture(tile);
    CloseWindow();    
    return EXIT_SUCCESS;
}

 
void show_title_screen(std::string_view title, std::string_view subtitle)
{
    int title_size       = MeasureText(title.data(), TITLE_FONT_SIZE);
    int subtitle_size    = MeasureText(subtitle.data(), SUBTITLE_FONT_SIZE);

    BeginDrawing();
        ClearBackground(BLACK);
        DrawText(title.data(),    (GetScreenWidth() - title_size   )/2, GetScreenHeight()*0.4,    TITLE_FONT_SIZE, WHITE);
        DrawText(subtitle.data(), (GetScreenWidth() - subtitle_size)/2, GetScreenHeight()*0.5, SUBTITLE_FONT_SIZE, WHITE);
    EndDrawing(); 
}
