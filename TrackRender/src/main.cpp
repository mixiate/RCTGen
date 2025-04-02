#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif
#include <math.h>

#include "track.h"

context_t get_context(light_t* lights,uint32_t num_lights,uint32_t dither)
{
	context_t context;
	context_init(&context,lights,num_lights,dither,palette_rct2(),TILE_SIZE);
	//context.palette.colors[0].r=0;
	//context.palette.colors[0].g=0;
	//context.palette.colors[0].b=0;
	return context;
}

int load_model(mesh_t* model,json_t* json,const char* name)
{
	json_t* mesh=json_object_get(json,name);
	if(mesh !=NULL)
	{
		if(json_is_string(mesh))
		{
			if(mesh_load_transform(model,json_string_value(mesh),rotate_y(-0.5*M_PI)))
			{
				printf("Failed to load model from file \"%s\"\n",json_string_value(mesh));
				return 1;
			}
			return 0;
		}
		printf("Error: Property \"%s\" not found or is not an object\n",name);
		return 1;
	}
	return 2;
}

int load_groups(json_t* json,uint64_t* out)
{
	//Load track sections
	uint64_t groups=0;
	for(int i=0; i<json_array_size(json); i++)
	{
		json_t* group_name=json_array_get(json,i);
		assert(group_name !=NULL);
		if(!json_is_string(group_name))
		{
			printf("Error: Array \"sections\" contains non-string value\n");
			return 1;
		}
		if(strcmp(json_string_value(group_name),"flat") ==0)groups|=TRACK_GROUP_FLAT;
		else if(strcmp(json_string_value(group_name),"brakes") ==0)groups|=TRACK_GROUP_BRAKES;
		else if(strcmp(json_string_value(group_name),"block_brakes") ==0)groups|=TRACK_GROUP_BLOCK_BRAKES;
		else if(strcmp(json_string_value(group_name),"diagonal_brakes") ==0)groups|=TRACK_GROUP_DIAGONAL_BRAKES;
		else if(strcmp(json_string_value(group_name),"sloped_brakes") ==0)groups|=TRACK_GROUP_SLOPED_BRAKES;
		else if(strcmp(json_string_value(group_name),"magnetic_brakes") ==0)groups|=TRACK_GROUP_MAGNETIC_BRAKES;
		else if(strcmp(json_string_value(group_name),"turns") ==0)groups|=TRACK_GROUP_TURNS;
		else if(strcmp(json_string_value(group_name),"gentle_slopes") ==0)groups|=TRACK_GROUP_GENTLE_SLOPES;
		else if(strcmp(json_string_value(group_name),"steep_slopes") ==0)groups|=TRACK_GROUP_STEEP_SLOPES;
		else if(strcmp(json_string_value(group_name),"vertical_slopes") ==0)groups|=TRACK_GROUP_VERTICAL_SLOPES;
		else if(strcmp(json_string_value(group_name),"diagonals") ==0)groups|=TRACK_GROUP_DIAGONALS;
		else if(strcmp(json_string_value(group_name),"sloped_turns") ==0)groups|=TRACK_GROUP_SLOPED_TURNS|TRACK_GROUP_STEEP_SLOPED_TURNS;
		else if(strcmp(json_string_value(group_name),"gentle_sloped_turns") ==0)groups|=TRACK_GROUP_SLOPED_TURNS;
		else if(strcmp(json_string_value(group_name),"banked_turns") ==0)groups|=TRACK_GROUP_BANKED_TURNS;
		else if(strcmp(json_string_value(group_name),"banked_sloped_turns") ==0)groups|=TRACK_GROUP_BANKED_SLOPED_TURNS;
		else if(strcmp(json_string_value(group_name),"large_sloped_turns") ==0)groups|=TRACK_GROUP_LARGE_SLOPED_TURNS;
		else if(strcmp(json_string_value(group_name),"large_banked_sloped_turns") ==0)groups|=TRACK_GROUP_LARGE_BANKED_SLOPED_TURNS;
		else if(strcmp(json_string_value(group_name),"s_bends") ==0)groups|=TRACK_GROUP_S_BENDS;
		else if(strcmp(json_string_value(group_name),"banked_s_bends") ==0)groups|=TRACK_GROUP_BANKED_S_BENDS;
		else if(strcmp(json_string_value(group_name),"helices") ==0)groups|=TRACK_GROUP_HELICES;
		else if(strcmp(json_string_value(group_name),"small_slope_transitions") ==0)groups|=TRACK_GROUP_SMALL_SLOPE_TRANSITIONS;
		else if(strcmp(json_string_value(group_name),"large_slope_transitions") ==0)groups|=TRACK_GROUP_LARGE_SLOPE_TRANSITIONS;
		else if(strcmp(json_string_value(group_name),"barrel_rolls") ==0)groups|=TRACK_GROUP_BARREL_ROLLS;
		else if(strcmp(json_string_value(group_name),"inline_twists") ==0)groups|=TRACK_GROUP_INLINE_TWISTS;
		else if(strcmp(json_string_value(group_name),"quarter_loops") ==0)groups|=TRACK_GROUP_QUARTER_LOOPS;
		else if(strcmp(json_string_value(group_name),"corkscrews") ==0)groups|=TRACK_GROUP_CORKSCREWS;
		else if(strcmp(json_string_value(group_name),"large_corkscrews") ==0)groups|=TRACK_GROUP_LARGE_CORKSCREWS;
		else if(strcmp(json_string_value(group_name),"half_loops") ==0)groups|=TRACK_GROUP_HALF_LOOPS;
		else if(strcmp(json_string_value(group_name),"vertical_loops")==0)groups|=TRACK_GROUP_VERTICAL_LOOPS;
		else if(strcmp(json_string_value(group_name),"medium_half_loops") ==0)groups|=TRACK_GROUP_MEDIUM_HALF_LOOPS;
		else if(strcmp(json_string_value(group_name),"large_half_loops") ==0)groups|=TRACK_GROUP_LARGE_HALF_LOOPS;
		else if(strcmp(json_string_value(group_name),"zero_g_rolls") ==0)groups|=TRACK_GROUP_ZERO_G_ROLLS;
		else if(strcmp(json_string_value(group_name),"dive_loops") ==0)groups|=TRACK_GROUP_DIVE_LOOPS;
		else if(strcmp(json_string_value(group_name),"boosters") ==0)groups|=TRACK_GROUP_BOOSTERS;
		else if(strcmp(json_string_value(group_name),"launched_lifts") ==0)groups|=TRACK_GROUP_LAUNCHED_LIFTS;
		else if(strcmp(json_string_value(group_name),"turn_bank_transitions") ==0)groups|=TRACK_GROUP_TURN_BANK_TRANSITIONS;
		else if(strcmp(json_string_value(group_name),"vertical_boosters") ==0)groups|=TRACK_GROUP_VERTICAL_BOOSTERS;
		else
		{
			printf("Error: Unrecognized section group \"%s\"\n",json_string_value(group_name));
			return 1;
		}
	}
	*out=groups;
	return 0;
}

int load_offsets(json_t* json,float* offsets)
{
const char* row_names[10]={"flat","gentle","steep","flat_banked","gentle_banked","inverted","diagonal","diagonal_banked","diagonal_gentle","diagonal_steep"};

//Zero offset array
memset(offsets,0,88*sizeof(float));

//Load offsets
	for(int i=0;i<10;i++)
	{
	json_t* row=json_object_get(json,row_names[i]);
		if(row == NULL)continue;
		if(!json_is_array(row) || json_array_size(row) != 8)
		{
		printf("Property \"%s\" is not an array of length 8\n",row_names[i]);
		return 1;
		}
		for(int j=0;j<8;j++)
		{
		json_t* value=json_array_get(row,j);
			if(!json_is_number(value))
			{
			printf("Array \"%s\" contains non numeric value\n",row_names[i]);
			return 1;
			}
		offsets[8*i+j]=json_number_value(value);
		}
	}
return 0;
}

int load_track_type(track_type_t* track_type,json_t* json)
{
	//Load track flags
	track_type->flags=0;
	json_t* flags=json_object_get(json,"flags");
	if(flags !=NULL)
	{
		if(!json_is_array(flags))
		{
			printf("Error: Property \"flags\" is not an array\n");
			return 1;
		}
		for(int i=0; i<json_array_size(flags); i++)
		{
			json_t* flag_name=json_array_get(flags,i);
			assert(flag_name !=NULL);
			if(!json_is_string(flag_name))
			{
				printf("Error: Array \"flags\" contains non-string value\n");
				return 1;
			}
			if(strcmp(json_string_value(flag_name),"has_lift") ==0)track_type->flags|=TRACK_HAS_LIFT;
			else if(strcmp(json_string_value(flag_name),"has_supports") ==0)track_type->flags|=TRACK_HAS_SUPPORTS;
			else if(strcmp(json_string_value(flag_name),"semi_split") ==0)track_type->flags|=TRACK_SEMI_SPLIT;
			else if(strcmp(json_string_value(flag_name),"split") ==0)track_type->flags|=TRACK_SPLIT;
			else if(strcmp(json_string_value(flag_name),"no_lift_sprite") ==0)track_type->flags|=TRACK_NO_LIFT_SPRITE;
			else if(strcmp(json_string_value(flag_name),"separate_tie") ==0)track_type->flags|=TRACK_SEPARATE_TIE;
			else if(strcmp(json_string_value(flag_name),"tie_at_boundary") ==0)track_type->flags|=TRACK_SEPARATE_TIE|TRACK_TIE_AT_BOUNDARY;
			else if(strcmp(json_string_value(flag_name),"special_end_offsets") ==0)track_type->flags|=TRACK_SPECIAL_OFFSETS;
			else
			{
				printf("Error: Unrecognized flag \"%s\"\n",json_string_value(flag_name));
				return 1;
			}
		}
	}

	json_t* groups=json_object_get(json,"sections");
	if(groups !=NULL)
	{
		if(!json_is_array(groups))
		{
			printf("Error: Property \"sections\" is not an array\n");
			return 1;
		}
		// Error reporting contained in load_groups
		if(load_groups(groups,&(track_type->groups)))return 1;
	}

	if(track_type->flags&TRACK_HAS_LIFT)
	{
		json_t* groups=json_object_get(json,"lift_sections");
		if(groups !=NULL)
		{
			if(!json_is_array(groups))
			{
				printf("Error: Property \"lift_sections\" is not an array\n");
				return 1;
			}
			if(load_groups(groups,&(track_type->lift_groups)))return 1;
		}
		json_t* offset=json_object_get(json,"lift_offset");
		if(offset)
		{
			if(!json_is_integer(offset))
			{
				printf("Error: Property \"lift_offset\" is not an int\n");
				return 1;
			}
			track_type->lift_offset=json_integer_value(offset);
		}
		else track_type->lift_offset=13;
	}

	//Load length
	json_t* length=json_object_get(json,"length");
	if(length !=NULL&&json_is_number(length))track_type->length=json_number_value(length)*TILE_SIZE;
	else
	{
		printf("Error: Property \"length\" not found or is not a number\n");
		return 1;
	}

	//Load brake length
	json_t* brake_length=json_object_get(json,"brake_length");
	if(brake_length !=NULL)
	{
		if(json_is_number(brake_length))track_type->brake_length=json_number_value(brake_length)*TILE_SIZE;
		else
		{
			printf("Error: Property \"brake_length\" not found or is not a number\n");
			return 1;
		}
	}
	else track_type->brake_length=TILE_SIZE;

	//Load tie length
	if(track_type->flags&TRACK_TIE_AT_BOUNDARY)
	{
		json_t* tie_length=json_object_get(json,"tie_length");
		if(tie_length !=NULL&&json_is_number(tie_length))track_type->tie_length=json_number_value(tie_length)*TILE_SIZE;
		else
		{
			printf("Error: Property \"tie_length\" not found or is not a number\n");
			return 1;
		}
	}

	//Load Z offset
	json_t* z_offset=json_object_get(json,"z_offset");
	if(z_offset !=NULL&&json_is_number(z_offset))track_type->z_offset=json_number_value(z_offset);
	else
	{
		printf("Error: Property \"z_offset\" not found or is not a number\n");
		return 1;
	}

	//Load support_spacing
	json_t* support_spacing=json_object_get(json,"support_spacing");
	if(support_spacing !=NULL)
	{
		if(json_is_number(support_spacing))track_type->support_spacing=json_number_value(support_spacing)*TILE_SIZE;
		else
		{
			printf("Error: Property \"support_spacing\" not found or is not a number\n");
			return 1;
		}
	}
	else track_type->support_spacing=TILE_SIZE;

	//Load pivot
	json_t* pivot=json_object_get(json,"pivot");
	if(pivot !=NULL)
	{
		if(json_is_number(pivot))track_type->pivot=json_number_value(pivot)*TILE_SIZE;
		else
		{
			printf("Error: Property \"pivot\" not found or is not a number\n");
			return 1;
		}
	}
	else track_type->pivot=0;

	//Load offset table
	if(track_type->flags & TRACK_SPECIAL_OFFSETS)
	{
	json_t* offsets=json_object_get(json,"offsets");
		if(offsets ==NULL || !json_is_object(offsets))
		{
		printf("Error: Property \"offsets\" not found or is not an object\n");
		return 1;
		}
		if(load_offsets(offsets,track_type->offset_table))return 1;
	}

	//Load models
	json_t* models=json_object_get(json,"models");
	if(models ==NULL||!json_is_object(models))
	{
		printf("Error: Property \"models\" not found or is not an object\n");
		return 1;
	}

	if(load_model(&(track_type->mesh),models,"track"))
	{
		printf("Error: Track mesh not found\n");
		return 1;
	}
	if(load_model(&(track_type->mask),models,"mask"))
	{
		mesh_destroy(&(track_type->mesh));
		printf("Error: Mask mesh not found\n");
		return 1;
	}

	if(track_type->flags&TRACK_HAS_LIFT)
	{
		if(load_model(&(track_type->lift_mesh),models,"lift"))
		{
			mesh_destroy(&(track_type->mesh));
			mesh_destroy(&(track_type->mask));
			printf("Error: Lift mesh not found\n");
			return 1;
		}
	}

	if(track_type->flags&TRACK_SEPARATE_TIE)
	{
		if(load_model(&(track_type->tie_mesh),models,"tie"))
		{
			mesh_destroy(&(track_type->mesh));
			mesh_destroy(&(track_type->mask));
			if(track_type->flags&TRACK_HAS_LIFT)mesh_destroy(&(track_type->lift_mesh));
			printf("Error: separate tie mesh not found\n");
			return 1;
		}

		if(track_type->flags&TRACK_TIE_AT_BOUNDARY)
		{
			if(load_model(&(track_type->mesh_tie),models,"track_tie"))
			{
				mesh_destroy(&(track_type->mesh));
				mesh_destroy(&(track_type->mask));
				mesh_destroy(&(track_type->tie_mesh));
				if(track_type->flags&TRACK_HAS_LIFT)mesh_destroy(&(track_type->lift_mesh));
				printf("Error: track_tie mesh not found\n");
				return 1;
			}
			if(track_type->flags&TRACK_HAS_LIFT)
			{
				if(load_model(&(track_type->lift_mesh_tie),models,"lift_tie"))
				{
					mesh_destroy(&(track_type->mesh));
					mesh_destroy(&(track_type->mask));
					mesh_destroy(&(track_type->tie_mesh));
					mesh_destroy(&(track_type->mesh_tie));
					mesh_destroy(&(track_type->lift_mesh));
					printf("Error: lift_tie mesh not found\n");
					return 1;
				}
			}
		}


	}

	const char* support_model_names[NUM_MODELS]={
	    "track_alt",
	    "support_flat",
	    "support_bank_sixth",
	    "support_bank_third",
	    "support_bank_half",
	    "support_bank_two_thirds",
	    "support_bank_five_sixths",
	    "support_bank",
	    "support_base",
	    "brake",
	    "block_brake",
	    "booster",
	    "magnetic_brake",
	    "support_steep_to_vertical",
	    "support_vertical_to_steep",
	    "support_vertical",
	    "support_vertical_twist",
	    "support_barrel_roll",
	    "support_half_loop",
	    "support_quarter_loop",
	    "support_corkscrew",
	    "support_zero_g_roll",
	    "support_large_zero_g_roll"};

	track_type->models_loaded=0;
	for(int i=0; i<NUM_MODELS; i++)
	{
		int result=load_model(&(track_type->models[i]),models,support_model_names[i]);
		if(result ==0)track_type->models_loaded|=1<<i;
		else if(result ==1)
		{
			mesh_destroy(&(track_type->mesh));
			mesh_destroy(&(track_type->mask));
			if(track_type->flags&TRACK_HAS_LIFT)mesh_destroy(&(track_type->lift_mesh));
			for(int j=0; j<i; j++)mesh_destroy(&(track_type->models[j]));
			printf("Error: failed to load model %s\n",support_model_names[i]);
			return 1;
		}
	}

	return 0;
}

int load_vector(vector3_t* vector,json_t* array)
{
	int size=json_array_size(array);
	if(size !=3)
	{
		printf("Vector must have 3 components\n");
		return 1;
	}

	json_t* x=json_array_get(array,0);
	json_t* y=json_array_get(array,1);
	json_t* z=json_array_get(array,2);

	if(!json_is_number(x)||!json_is_number(y)||!json_is_number(z))
	{
		printf("Vector components must be numeric\n");
		return 1;
	}
	vector->x=json_number_value(x);
	vector->y=json_number_value(y);
	vector->z=json_number_value(z);
	return 0;
}

int load_lights(light_t* lights,int* lights_count,json_t* json)
{
	int num_lights=json_array_size(json);
	for(int i=0; i<num_lights; i++)
	{
		json_t* light=json_array_get(json,i);
		assert(light !=NULL);
		if(!json_is_object(light))
		{
			printf("Warning: Light array contains an element which is not an object-ignoring\n");
			continue;
		}

		json_t* type=json_object_get(light,"type");
		if(type ==NULL||!json_is_string(type))
		{
			printf("Error: Property \"type\" not found or is not a string\n");
			return 1;
		}

		const char* type_value=json_string_value(type);
		if(strcmp(type_value,"diffuse") ==0)lights[i].type=LIGHT_DIFFUSE;
		else if(strcmp(type_value,"specular") ==0)lights[i].type=LIGHT_SPECULAR;
		else
		{
			printf("Unrecognized light type \"%s\"\n",type);
			free(lights);
		}

		json_t* shadow=json_object_get(light,"shadow");
		if(shadow ==NULL||!json_is_boolean(shadow))
		{
			printf("Error: Property \"shadow\" not found or is not a boolean\n");
			return 1;
		}
		if(json_boolean_value(shadow))lights[i].shadow=1;
		else lights[i].shadow=0;

		json_t* direction=json_object_get(light,"direction");
		if(direction ==NULL||!json_is_array(direction))
		{
			printf("Error: Property \"direction\" not found or is not a direction\n");
			return 1;
		}
		if(load_vector(&(lights[i].direction),direction))return 1;
		lights[i].direction=vector3_normalize(lights[i].direction);

		json_t* strength=json_object_get(light,"strength");
		if(strength ==NULL||!json_is_number(strength))
		{
			printf("Error: Property \"strength\" not found or is not a number\n");
			return 1;
		}
		lights[i].intensity=json_number_value(strength);
	}
	*lights_count=num_lights;
	return 0;
}

int is_in_mask(int x, int y, mask_t* mask);

#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
#include <format>

void dump_mask(const char* name, track_section_t& track_section) {
	char desc_file_name[512];
	sprintf(desc_file_name, "C:/Files/masks/default/%s.json", name);
	std::ofstream mask_desc;
	mask_desc.open(desc_file_name);

	mask_desc << "[\n";

	for (int view_i = 0; view_i < 4; view_i++) {
		const auto view = track_section.views[view_i];

		const int dimensions = 768;

		std::vector<std::tuple<std::string, image_t>> images;

		mask_desc << "    {\n";

		mask_desc << std::format("        \"uses_track_mask\": {},\n", view.flags & VIEW_NEEDS_TRACK_MASK ? true : false);
		mask_desc << std::format("        \"enforce_non_overlapping\": {},\n", view.flags & VIEW_ENFORCE_NON_OVERLAPPING ? true : false);
		mask_desc << "        \"masks\": [\n";

		if (view.num_sprites == 0)
		{
			mask_desc << "        ]\n";

			if (view_i == 3) {
				mask_desc << "    }\n";
			}
			else {
				mask_desc << "    },\n";
			}

			continue;
		}

		for (int mask_i = 0; mask_i < view.num_sprites; mask_i++)
		{
			mask_desc << "            {\n";

			if (view.masks == nullptr) {
				mask_desc << "                \"use_mask\": false,\n";
				mask_desc << "                \"image_path\": \"\",\n";
				mask_desc << "                \"color\": 0,\n";
				mask_desc << "                \"mask_type\": \"none\",\n";
				mask_desc << "                \"offset_x\": " << 0 << ",\n";
				mask_desc << "                \"offset_y\": " << 0 << ",\n";
				mask_desc << "                \"flipped\": false\n";
			}
			else {
				auto mask = view.masks[mask_i];

				image_t image;
				image_new(&image, dimensions, dimensions, 0, 0, 0);

				for (int x = 0; x < dimensions; x++) {
					for (int y = 0; y < dimensions; y++) {
						image.pixels[y * dimensions + x] = is_in_mask(x - (dimensions / 2), y - (dimensions / 2), &mask) ? 202 : 0;
					}
				}

				bool unique = false;

				if (images.size() > 0) {
					const auto other_image = std::get<1>(images.at(images.size() - 1));

					size_t first_sum = 0;
					size_t second_sum = 0;

					for (int x = 0; x < dimensions; x++) {
						for (int y = 0; y < dimensions; y++) {
							first_sum += image.pixels[y * dimensions + x];
							second_sum += other_image.pixels[y * dimensions + x];

							if (first_sum != second_sum) {
								unique = true;
							}
						}
					}
				}

				int color = 202 + images.size();
				if (images.size() && !unique) {
					color -= 1;
				}

				mask_desc << "                \"use_mask\": true,\n";

				char image_path_buffer[512];
				sprintf(image_path_buffer, "                \"image_path\": \"%s_%d.png\",\n", name, view_i + 1);

				mask_desc << image_path_buffer;
				mask_desc << std::format("                \"color\": {},\n", color);
				if (mask.track_mask_op == TRACK_MASK_NONE) {
					mask_desc << "                \"mask_type\": \"none\",\n";
				}
				if (mask.track_mask_op == TRACK_MASK_DIFFERENCE) {
					mask_desc << "                \"mask_type\": \"difference\",\n";
				}
				if (mask.track_mask_op == TRACK_MASK_INTERSECT) {
					mask_desc << "                \"mask_type\": \"intersect\",\n";
				}
				if (mask.track_mask_op == TRACK_MASK_UNION) {
					mask_desc << "                \"mask_type\": \"union\",\n";
				}
				if (mask.track_mask_op == TRACK_MASK_TRANSFER_NEXT) {
					mask_desc << "                \"mask_type\": \"transfer_next\",\n";
				}
				mask_desc << "                \"offset_x\": " << mask.x_offset << ",\n";
				mask_desc << "                \"offset_y\": " << mask.y_offset << ",\n";
				mask_desc << "                \"flipped\": false\n";

				if (unique || images.size() == 0) {
					const auto image_file_path = std::format("C:/Files/masks/default/{}_{}_{}.png", name, view_i + 1, images.size() + 1);

					images.push_back(std::make_tuple(image_file_path, image));
				}
			}

			if (mask_i == view.num_sprites - 1) {
				mask_desc << "            }\n";
			}
			else {
				mask_desc << "            },\n";
			}
		}

		mask_desc << "        ]\n";

		if (view_i == 3) {
			mask_desc << "    }\n";
		}
		else {
			mask_desc << "    },\n";
		}

		bool overlaps = false;

		for (size_t i = 0; i < images.size(); i++) {
			const image_t& current_image = std::get<1>(images.at(i));

			for (size_t j = 0; j < images.size(); j++) {
				if (i == j)
					continue;

				const image_t& other_image = std::get<1>(images.at(j));

				for (int x = 0; x < dimensions; x++) {
					for (int y = 0; y < dimensions; y++) {
						if (current_image.pixels[y * dimensions + x] != 0 && other_image.pixels[y * dimensions + x] != 0) {
							overlaps = true;
							goto break_out;
						}
					}
				}
			}
		}

	break_out:
		if (images.size()) {
			if (overlaps) {
				printf("%s\n", name);
				for (auto image_desc : images) {
					FILE* file;
					file = fopen(std::get<0>(image_desc).c_str(), "wb");

					image_write_png(&std::get<1>(image_desc), file);

					fclose(file);
				}
			}
			else {
				image_t image;
				image_new(&image, dimensions, dimensions, 0, 0, 0);

				for (size_t i = 0; i < images.size(); i++) {
					const image_t& split_image = std::get<1>(images.at(i));

					for (int x = 0; x < dimensions; x++) {
						for (int y = 0; y < dimensions; y++) {
							if (split_image.pixels[y * dimensions + x]) {
								image.pixels[y * dimensions + x] = split_image.pixels[y * dimensions + x] + i;
							}
						}
					}
				}

				const auto image_file_path = std::format("C:/Files/masks/default/{}_{}.png", name, view_i + 1);

				FILE* file;
				file = fopen(image_file_path.c_str(), "wb");

				image_write_png(&image, file);

				fclose(file);
			}
		}
	}

	mask_desc << "]\n";

	mask_desc.close();
}

void dump_masks() {
	track_list_t track_list = track_list_default;

	/*dump_mask("flat", track_list.flat);
	dump_mask("flat_asymmetric", track_list.flat_asymmetric);
	dump_mask("brake", track_list.brake);
	dump_mask("brake_diag", track_list.brake_diag);
	dump_mask("brake_gentle", track_list.brake_gentle);
	dump_mask("brake_gentle_diag", track_list.brake_gentle_diag);
	dump_mask("magnetic_brake", track_list.magnetic_brake);
	dump_mask("magnetic_brake_diag", track_list.magnetic_brake_diag);
	dump_mask("magnetic_brake_gentle", track_list.magnetic_brake_gentle);
	dump_mask("magnetic_brake_gentle_diag", track_list.magnetic_brake_gentle_diag);
	dump_mask("block_brake", track_list.block_brake);
	dump_mask("block_brake_diag", track_list.block_brake_diag);
	dump_mask("booster", track_list.booster);
	dump_mask("flat_to_gentle_up", track_list.flat_to_gentle_up);
	dump_mask("gentle_up_to_flat", track_list.gentle_up_to_flat);
	dump_mask("gentle", track_list.gentle);
	dump_mask("gentle_to_steep_up", track_list.gentle_to_steep_up);
	dump_mask("steep_to_gentle_up", track_list.steep_to_gentle_up);
	dump_mask("steep", track_list.steep);
	dump_mask("steep_to_vertical_up", track_list.steep_to_vertical_up);
	dump_mask("vertical_to_steep_up", track_list.vertical_to_steep_up);
	dump_mask("vertical", track_list.vertical);
	dump_mask("small_turn_left", track_list.small_turn_left);
	dump_mask("medium_turn_left", track_list.medium_turn_left);
	dump_mask("large_turn_left_to_diag", track_list.large_turn_left_to_diag);
	dump_mask("large_turn_right_to_diag", track_list.large_turn_right_to_diag);
	dump_mask("flat_diag", track_list.flat_diag);
	dump_mask("flat_to_gentle_up_diag", track_list.flat_to_gentle_up_diag);
	dump_mask("gentle_to_flat_up_diag", track_list.gentle_to_flat_up_diag);
	dump_mask("gentle_diag", track_list.gentle_diag);
	dump_mask("gentle_to_steep_up_diag", track_list.gentle_to_steep_up_diag);
	dump_mask("steep_to_gentle_up_diag", track_list.steep_to_gentle_up_diag);
	dump_mask("steep_diag", track_list.steep_diag);
	dump_mask("flat_to_left_bank", track_list.flat_to_left_bank);
	dump_mask("flat_to_right_bank", track_list.flat_to_right_bank);
	dump_mask("left_bank_to_gentle_up", track_list.left_bank_to_gentle_up);
	dump_mask("right_bank_to_gentle_up", track_list.right_bank_to_gentle_up);
	dump_mask("gentle_up_to_left_bank", track_list.gentle_up_to_left_bank);
	dump_mask("gentle_up_to_right_bank", track_list.gentle_up_to_right_bank);
	dump_mask("left_bank", track_list.left_bank);
	dump_mask("flat_to_left_bank_diag", track_list.flat_to_left_bank_diag);
	dump_mask("flat_to_right_bank_diag", track_list.flat_to_right_bank_diag);
	dump_mask("left_bank_to_gentle_up_diag", track_list.left_bank_to_gentle_up_diag);
	dump_mask("right_bank_to_gentle_up_diag", track_list.right_bank_to_gentle_up_diag);
	dump_mask("gentle_up_to_left_bank_diag", track_list.gentle_up_to_left_bank_diag);
	dump_mask("gentle_up_to_right_bank_diag", track_list.gentle_up_to_right_bank_diag);
	dump_mask("left_bank_diag", track_list.left_bank_diag);
	dump_mask("small_turn_left_bank", track_list.small_turn_left_bank);
	dump_mask("medium_turn_left_bank", track_list.medium_turn_left_bank);
	dump_mask("large_turn_left_to_diag_bank", track_list.large_turn_left_to_diag_bank);
	dump_mask("large_turn_right_to_diag_bank", track_list.large_turn_right_to_diag_bank);
	dump_mask("small_turn_left_gentle_up", track_list.small_turn_left_gentle_up);
	dump_mask("small_turn_right_gentle_up", track_list.small_turn_right_gentle_up);
	dump_mask("medium_turn_left_gentle_up", track_list.medium_turn_left_gentle_up);
	dump_mask("medium_turn_right_gentle_up", track_list.medium_turn_right_gentle_up);
	dump_mask("very_small_turn_left_steep_up", track_list.very_small_turn_left_steep_up);
	dump_mask("very_small_turn_right_steep_up", track_list.very_small_turn_right_steep_up);
	dump_mask("vertical_twist_left_up", track_list.vertical_twist_left_up);
	dump_mask("vertical_twist_right_up", track_list.vertical_twist_right_up);
	dump_mask("gentle_up_to_gentle_up_left_bank", track_list.gentle_up_to_gentle_up_left_bank);
	dump_mask("gentle_up_to_gentle_up_right_bank", track_list.gentle_up_to_gentle_up_right_bank);
	dump_mask("gentle_up_left_bank_to_gentle_up", track_list.gentle_up_left_bank_to_gentle_up);
	dump_mask("gentle_up_right_bank_to_gentle_up", track_list.gentle_up_right_bank_to_gentle_up);
	dump_mask("left_bank_to_gentle_up_left_bank", track_list.left_bank_to_gentle_up_left_bank);
	dump_mask("gentle_up_left_bank_to_left_bank", track_list.gentle_up_left_bank_to_left_bank);
	dump_mask("right_bank_to_gentle_up_right_bank", track_list.right_bank_to_gentle_up_right_bank);
	dump_mask("gentle_up_right_bank_to_right_bank", track_list.gentle_up_right_bank_to_right_bank);
	dump_mask("gentle_up_left_bank", track_list.gentle_up_left_bank);
	dump_mask("gentle_up_right_bank", track_list.gentle_up_right_bank);
	dump_mask("flat_to_gentle_up_left_bank", track_list.flat_to_gentle_up_left_bank);
	dump_mask("flat_to_gentle_up_right_bank", track_list.flat_to_gentle_up_right_bank);
	dump_mask("gentle_up_left_bank_to_flat", track_list.gentle_up_left_bank_to_flat);
	dump_mask("gentle_up_right_bank_to_flat", track_list.gentle_up_right_bank_to_flat);
	dump_mask("small_turn_left_bank_gentle_up", track_list.small_turn_left_bank_gentle_up);
	dump_mask("small_turn_right_bank_gentle_up", track_list.small_turn_right_bank_gentle_up);
	dump_mask("medium_turn_left_bank_gentle_up", track_list.medium_turn_left_bank_gentle_up);
	dump_mask("medium_turn_right_bank_gentle_up", track_list.medium_turn_right_bank_gentle_up);
	dump_mask("s_bend_left", track_list.s_bend_left);
	dump_mask("s_bend_right", track_list.s_bend_right);
	dump_mask("s_bend_left_bank", track_list.s_bend_left_bank);
	dump_mask("s_bend_right_bank", track_list.s_bend_right_bank);
	dump_mask("small_helix_left_up", track_list.small_helix_left_up);
	dump_mask("small_helix_right_up", track_list.small_helix_right_up);
	dump_mask("medium_helix_left_up", track_list.medium_helix_left_up);
	dump_mask("medium_helix_right_up", track_list.medium_helix_right_up);
	dump_mask("barrel_roll_left", track_list.barrel_roll_left);
	dump_mask("barrel_roll_right", track_list.barrel_roll_right);
	dump_mask("inline_twist_left", track_list.inline_twist_left);
	dump_mask("inline_twist_right", track_list.inline_twist_right);
	dump_mask("half_loop", track_list.half_loop);
	dump_mask("vertical_loop_left", track_list.left_vertical_loop);
	dump_mask("vertical_loop_right", track_list.right_vertical_loop);
	dump_mask("medium_half_loop_left", track_list.medium_half_loop_left);
	dump_mask("medium_half_loop_right", track_list.medium_half_loop_right);
	dump_mask("large_half_loop_left", track_list.large_half_loop_left);
	dump_mask("large_half_loop_right", track_list.large_half_loop_right);
	dump_mask("flat_to_steep_up", track_list.flat_to_steep_up);
	dump_mask("steep_to_flat_up", track_list.steep_to_flat_up);
	dump_mask("small_flat_to_steep_up", track_list.small_flat_to_steep_up);
	dump_mask("small_steep_to_flat_up", track_list.small_steep_to_flat_up);
	dump_mask("small_flat_to_steep_up_diag", track_list.small_flat_to_steep_up_diag);
	dump_mask("small_steep_to_flat_up_diag", track_list.small_steep_to_flat_up_diag);
	dump_mask("quarter_loop_up", track_list.quarter_loop_up);
	dump_mask("corkscrew_left", track_list.corkscrew_left);
	dump_mask("corkscrew_right", track_list.corkscrew_right);
	dump_mask("large_corkscrew_left", track_list.large_corkscrew_left);
	dump_mask("large_corkscrew_right", track_list.large_corkscrew_right);
	dump_mask("zero_g_roll_left", track_list.zero_g_roll_left);
	dump_mask("zero_g_roll_right", track_list.zero_g_roll_right);
	dump_mask("large_zero_g_roll_left", track_list.large_zero_g_roll_left);
	dump_mask("large_zero_g_roll_right", track_list.large_zero_g_roll_right);
	dump_mask("small_turn_left_bank_to_gentle_up", track_list.small_turn_left_bank_to_gentle_up);
	dump_mask("small_turn_right_bank_to_gentle_up", track_list.small_turn_right_bank_to_gentle_up);
	dump_mask("launched_lift", track_list.launched_lift);
	dump_mask("large_turn_left_to_diag_gentle_up", track_list.large_turn_left_to_diag_gentle_up);
	dump_mask("large_turn_right_to_diag_gentle_up", track_list.large_turn_right_to_diag_gentle_up);
	dump_mask("large_turn_left_to_orthogonal_gentle_up", track_list.large_turn_left_to_orthogonal_gentle_up);
	dump_mask("large_turn_right_to_orthogonal_gentle_up", track_list.large_turn_right_to_orthogonal_gentle_up);
	dump_mask("gentle_up_to_gentle_up_left_bank_diag", track_list.gentle_up_to_gentle_up_left_bank_diag);
	dump_mask("gentle_up_to_gentle_up_right_bank_diag", track_list.gentle_up_to_gentle_up_right_bank_diag);
	dump_mask("gentle_up_left_bank_to_gentle_up_diag", track_list.gentle_up_left_bank_to_gentle_up_diag);
	dump_mask("gentle_up_right_bank_to_gentle_up_diag", track_list.gentle_up_right_bank_to_gentle_up_diag);
	dump_mask("left_bank_to_gentle_up_left_bank_diag", track_list.left_bank_to_gentle_up_left_bank_diag);
	dump_mask("right_bank_to_gentle_up_right_bank_diag", track_list.right_bank_to_gentle_up_right_bank_diag);
	dump_mask("gentle_up_left_bank_to_left_bank_diag", track_list.gentle_up_left_bank_to_left_bank_diag);
	dump_mask("gentle_up_right_bank_to_right_bank_diag", track_list.gentle_up_right_bank_to_right_bank_diag);
	dump_mask("gentle_up_left_bank_diag", track_list.gentle_up_left_bank_diag);
	dump_mask("gentle_up_right_bank_diag", track_list.gentle_up_right_bank_diag);
	dump_mask("flat_to_gentle_up_left_bank_diag", track_list.flat_to_gentle_up_left_bank_diag);
	dump_mask("flat_to_gentle_up_right_bank_diag", track_list.flat_to_gentle_up_right_bank_diag);
	dump_mask("gentle_up_left_bank_to_flat_diag", track_list.gentle_up_left_bank_to_flat_diag);
	dump_mask("gentle_up_right_bank_to_flat_diag", track_list.gentle_up_right_bank_to_flat_diag);
	dump_mask("large_turn_left_bank_to_diag_gentle_up", track_list.large_turn_left_bank_to_diag_gentle_up);
	dump_mask("large_turn_right_bank_to_diag_gentle_up", track_list.large_turn_right_bank_to_diag_gentle_up);
	dump_mask("large_turn_left_bank_to_orthogonal_gentle_up", track_list.large_turn_left_bank_to_orthogonal_gentle_up);
	dump_mask("large_turn_right_bank_to_orthogonal_gentle_up", track_list.large_turn_right_bank_to_orthogonal_gentle_up);
	dump_mask("flat_to_steep_up_diag", track_list.flat_to_steep_up_diag);
	dump_mask("steep_to_flat_up_diag", track_list.steep_to_flat_up_diag);
	dump_mask("steep_to_vertical_up_diag", track_list.steep_to_vertical_up_diag);
	dump_mask("vertical_to_steep_up_diag", track_list.vertical_to_steep_up_diag);
	dump_mask("vertical_diag", track_list.vertical_diag);
	dump_mask("vertical_twist_left_to_diag_up", track_list.vertical_twist_left_to_diag_up);
	dump_mask("vertical_twist_right_to_diag_up", track_list.vertical_twist_right_to_diag_up);
	dump_mask("vertical_twist_left_to_orthogonal_up", track_list.vertical_twist_left_to_orthogonal_up);
	dump_mask("vertical_twist_right_to_orthogonal_up", track_list.vertical_twist_right_to_orthogonal_up);
	dump_mask("vertical_booster", track_list.vertical_booster);*/
	dump_mask("flat_to_steep_up_diag", track_list.flat_to_steep_up_diag);
	dump_mask("steep_to_flat_up_diag", track_list.steep_to_flat_up_diag);
	dump_mask("dive_loop_45_left", track_list.dive_loop_45_left);
	dump_mask("dive_loop_45_right", track_list.dive_loop_45_right);
	dump_mask("dive_loop_90_left", track_list.dive_loop_90_left);
	dump_mask("dive_loop_90_right", track_list.dive_loop_90_right);
}

int main(int argc,char** argv)
{
	dump_masks();
	return 0;
}
