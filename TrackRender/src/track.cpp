#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#define _CRT_SECURE_NO_WARNINGS
#endif
#include "track.h"
#include "sprites.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <array>
#include <format>

vector3_t start_offset;
vector3_t end_offset;

vector3_t change_coordinates(vector3_t a)
{
	vector3_t result={a.z,a.y,a.x};
	return result;
}

typedef struct
{
	track_point_t (*track_curve)(float);
	float scale;
	float offset;
	float y_offset;
	float z_offset;
	float length;
	int flags;
}track_transform_args_t;

track_point_t get_track_point(track_point_t (*curve)(float distance),int flags,float z_offset,float length,float u)
{
	track_point_t track_point;
	if(u<0)
	{
		track_point=curve(0);
		track_point.position=vector3_add(track_point.position,vector3_mult(track_point.tangent,u));
	}
	else if(u>length)
	{
		track_point=curve(length);
		track_point.position=vector3_add(track_point.position,vector3_mult(track_point.tangent,(u-length)));
	}
	else
	{
		track_point=curve(u);
	}

	if(flags&TRACK_DIAGONAL)track_point.position.x+=0.5*TILE_SIZE;
	if(flags&TRACK_DIAGONAL_2)track_point.position.z+=0.5*TILE_SIZE;
	track_point.position.y+=z_offset-2*CLEARANCE_HEIGHT;
	if(!(flags&TRACK_VERTICAL))track_point.position.z-=0.5*TILE_SIZE;

	float v=u/length;
	if(v<0)v=0;
	else if(v>1)v=1;

	track_point.position=vector3_add(track_point.position,vector3_mult(start_offset,2*v*v*v-3*v*v+1));
	track_point.position=vector3_add(track_point.position,vector3_mult(end_offset,-2*v*v*v+3*v*v));
	return track_point;
}

track_point_t only_yaw(track_point_t input)
{
	track_point_t output;
	output.position=input.position;
	output.normal=vector3(0,1,0);
	output.binormal=vector3_normalize(vector3_cross(output.normal,change_coordinates(input.tangent)));
	output.tangent=vector3_normalize(vector3_cross(output.normal,output.binormal));
	return output;
}

vertex_t track_transform(vector3_t vertex,vector3_t normal,void* data)
{
	track_transform_args_t args=*((track_transform_args_t*)data);

	vertex.z=args.scale*vertex.z+args.offset;

	track_point_t track_point=get_track_point(args.track_curve,args.flags,args.z_offset,args.length,vertex.z);
	if (args.y_offset != 0.0) 
		track_point.position = vector3_add(track_point.position, vector3_mult(track_point.binormal, args.y_offset));

	vertex_t out;
	out.vertex=change_coordinates(vector3_add(track_point.position,vector3_add(vector3_mult(track_point.normal,vertex.y),vector3_mult(track_point.binormal,vertex.x))));
	out.normal=change_coordinates(vector3_add(vector3_mult(track_point.tangent,normal.z),vector3_add(vector3_mult(track_point.normal,normal.y),vector3_mult(track_point.binormal,normal.x))));
	out.distance = vertex.z / args.length;
	return out;
}

vertex_t base_transform(vector3_t vertex,vector3_t normal,void* data)
{
	track_transform_args_t args=*((track_transform_args_t*)data);

	vertex.z=args.scale*vertex.z+args.offset;

	track_point_t track_point=get_track_point(args.track_curve,args.flags,args.z_offset,args.length,vertex.z);

	track_point.binormal=vector3_normalize(vector3_cross(vector3(0,1,0),track_point.tangent));
	track_point.normal=vector3_normalize(vector3_cross(track_point.tangent,track_point.binormal));

	vertex_t out;
	out.vertex=change_coordinates(vector3_add(track_point.position,vector3_add(vector3_mult(vector3(0,1,0),vertex.y),vector3_mult(track_point.binormal,vertex.x))));
	out.normal=change_coordinates(vector3_add(vector3_mult(track_point.tangent,normal.z),vector3_add(vector3_mult(track_point.normal,normal.y),vector3_mult(track_point.binormal,normal.x))));
	out.distance = vertex.z / args.length;
	return out;
}

#define DENOM 6

int get_support_index(int bank)
{
	switch(bank)
	{
	case 0:
		return MODEL_FLAT;
		break;
	case -1:
	case 1:
		return MODEL_BANK_SIXTH;
		break;
	case -2:
	case 2:
		return MODEL_BANK_THIRD;
	case -3:
	case 3:
		return MODEL_BANK_HALF;
		break;
	case -4:
	case 4:
		return MODEL_BANK_TWO_THIRDS;
		break;
	case -5:
	case 5:
		return MODEL_BANK_FIVE_SIXTHS;
	case -6:
	case 6:
		return MODEL_BANK;
		break;
	}
	return MODEL_FLAT;
}

int get_special_index(int flags)
{
	switch(flags&TRACK_SPECIAL_MASK)
	{
	case TRACK_SPECIAL_STEEP_TO_VERTICAL:
		return MODEL_SPECIAL_STEEP_TO_VERTICAL;
		break;
	case TRACK_SPECIAL_VERTICAL_TO_STEEP:
		return MODEL_SPECIAL_VERTICAL_TO_STEEP;
		break;
	case TRACK_SPECIAL_VERTICAL:
		return MODEL_SPECIAL_VERTICAL;
		break;
	case TRACK_SPECIAL_VERTICAL_TWIST_LEFT:
	case TRACK_SPECIAL_VERTICAL_TWIST_RIGHT:
		return MODEL_SPECIAL_VERTICAL_TWIST;
		break;
	case TRACK_SPECIAL_BARREL_ROLL_LEFT:
	case TRACK_SPECIAL_BARREL_ROLL_RIGHT:
		return MODEL_SPECIAL_BARREL_ROLL;
		break;
	case TRACK_SPECIAL_HALF_LOOP:
		return MODEL_SPECIAL_HALF_LOOP;
		break;
	case TRACK_SPECIAL_QUARTER_LOOP:
		return MODEL_SPECIAL_QUARTER_LOOP;
		break;
	case TRACK_SPECIAL_CORKSCREW_LEFT:
	case TRACK_SPECIAL_CORKSCREW_RIGHT:
		return MODEL_SPECIAL_CORKSCREW;
		break;
	case TRACK_SPECIAL_ZERO_G_ROLL_LEFT:
	case TRACK_SPECIAL_ZERO_G_ROLL_RIGHT:
		return MODEL_SPECIAL_ZERO_G_ROLL;
		break;
	case TRACK_SPECIAL_LARGE_ZERO_G_ROLL_LEFT:
	case TRACK_SPECIAL_LARGE_ZERO_G_ROLL_RIGHT:
		return MODEL_SPECIAL_LARGE_ZERO_G_ROLL;
		break;
	case TRACK_SPECIAL_BRAKE:
		return MODEL_SPECIAL_BRAKE;
		break;
	case TRACK_SPECIAL_BLOCK_BRAKE:
		return MODEL_SPECIAL_BLOCK_BRAKE;
		break;
	case TRACK_SPECIAL_MAGNETIC_BRAKE:
		return MODEL_SPECIAL_MAGNETIC_BRAKE;
		break;
	case TRACK_SPECIAL_BOOSTER:
	case TRACK_SPECIAL_LAUNCHED_LIFT:
	case TRACK_SPECIAL_VERTICAL_BOOSTER:
		return MODEL_SPECIAL_BOOSTER;
		break;
	}
	assert(0);
	return 0;
}

float recalculate_y_offset_track_section_length(track_section_t* track_section, float y_offset, float z_offset)
{
	if (y_offset == 0.0)
		return track_section->length;
		
	const size_t position_count = 5000;

	std::vector<vector3_t> track_section_positions(position_count);

	for (size_t i = 0; i < position_count; i++)
	{
		const track_point_t track_point = get_track_point(track_section->curve, track_section->flags, z_offset, track_section->length, (track_section->length / position_count) * i);
		auto offset_track_position = vector3_add(track_point.position, vector3_mult(track_point.binormal, y_offset));

		// reverse these calculations
		if (track_section->flags & TRACK_DIAGONAL)offset_track_position.x -= 0.5 * TILE_SIZE;
		if (track_section->flags & TRACK_DIAGONAL_2)offset_track_position.z -= 0.5 * TILE_SIZE;
		offset_track_position.y -= z_offset - 2 * CLEARANCE_HEIGHT;
		if (!(track_section->flags & TRACK_VERTICAL))offset_track_position.z += 0.5 * TILE_SIZE;

		track_section_positions.push_back(offset_track_position);
	}

	float distance = 0.0;
	vector3_t previous_position = track_section_positions.front();
	for (const auto position : track_section_positions) {
		const auto v = vector3_sub(previous_position, position);
		distance += std::sqrtf(vector3_dot(v, v));
		previous_position = position;
	}

	return distance - std::abs(y_offset);
}

void add_track_models(
	context_t* context,
	track_section_t* track_section,
	track_type_t* track_type,
	const int extrude_behind,
	const int extrude_in_front,
	const int track_mask,
	const int rendered_views,
	image_t* images,
	const int subtype,
	track_mesh_t& track_mesh)
{
	float track_section_length_end_addition = 0.0;

	for (auto& track_section_length_entry : track_type->track_section_length_entries)
	{
		if (track_section_length_entry.track_section_name == track_section->name)
		{
			track_section_length_end_addition = track_section_length_entry.end_addition;
			break;
		}
	}
	const float track_section_length = track_section->length + track_section_length_end_addition;

	float z_offset = ((track_type->z_offset / 8.0) * CLEARANCE_HEIGHT);

	float y_offset_track_section_length = recalculate_y_offset_track_section_length(track_section, track_mesh.y_offset, z_offset);

	int num_meshes = (int)floor(0.5 + y_offset_track_section_length / track_type->length);
	float scale = track_section_length / (num_meshes * track_type->length);

	//If alternating track meshes are used,we would prefer to have an even number of meshes as long as it doesn't cause too much distortion
	if (track_type->models_loaded & (1 << MODEL_TRACK_ALT))
	{
		int num_meshes_even = 2 * (int)floor(0.5 + track_section_length / (2 * track_type->length));
		if (track_section->flags & TRACK_ALT_PREFER_ODD)num_meshes_even = 2 * (int)floor(track_section_length / (2 * track_type->length)) + 1;
		float scale_even = track_section_length / (num_meshes_even * track_type->length);
		if (scale_even > 0.9 && scale_even < 1.11111)
		{
			num_meshes = num_meshes_even;
			scale = scale_even;
		}
	}

	float length = scale * track_type->length;

	mesh_t* mesh;
	mesh_t* mesh_tie;
	switch (subtype)
	{
	case TRACK_SUBTYPE_DEFAULT:
		mesh = &(track_mesh.mesh);
		mesh_tie = &(track_type->mesh_tie);
		break;
	case TRACK_SUBTYPE_LIFT:
		mesh = &(track_type->lift_mesh);
		mesh_tie = &(track_type->lift_mesh_tie);
		break;
	default:
		assert(0);
		break;
	}

	//Add ghost models/track masks at start and end
	const int mask_model_count = track_type->length < 0.5 ? 2 : 1;

	track_transform_args_t args;
	args.scale = scale;
	args.offset = -length;
	args.y_offset = 0.0;
	args.z_offset = z_offset;
	args.track_curve = track_section->curve;
	args.flags = track_section->flags;
	args.length = track_section_length;
	for (int i = 0; i < mask_model_count; i++)
	{
		args.offset = -length * float(i + 1);
		if (track_mask)
		{
			args.y_offset = 0.0;
			context_add_model_transformed(context, &(track_type->mask), track_transform, &args, 0, vector2(0.0, 0.0));
		}
		else if (!extrude_behind)
		{
			args.y_offset = track_mesh.y_offset;
			context_add_model_transformed(context, mesh, track_transform, &args, MESH_GHOST, vector2(0.0, 0.0));
		}

		args.offset = track_section_length + (length * float(i));
		if (track_mask)
		{
			args.y_offset = 0.0;
			context_add_model_transformed(context, &(track_type->mask), track_transform, &args, 0, vector2(0.0, 0.0));
		}
		else if (!extrude_in_front)
		{
			args.y_offset = track_mesh.y_offset;
			context_add_model_transformed(context, mesh, track_transform, &args, MESH_GHOST, vector2(0.0, 0.0));
		}
	}

	if (track_type->flags & TRACK_TIE_AT_BOUNDARY)
	{
		//Determine start angle
		int start_angle = 3;
		if (rendered_views & 1)start_angle = 0;
		else if (rendered_views & 2)start_angle = 1;
		else if (rendered_views & 4)start_angle = 2;
		if (track_section->flags & (TRACK_DIAGONAL | TRACK_DIAGONAL_2))start_angle = (start_angle + 1) % 4;
		int end_angle = start_angle;
		if (track_section->flags & TRACK_EXIT_90_DEG_LEFT)end_angle -= 1;
		else if (track_section->flags & TRACK_EXIT_90_DEG_RIGHT)end_angle += 1;
		else if (track_section->flags & TRACK_EXIT_180_DEG)end_angle += 2;
		else if ((track_section->flags & TRACK_EXIT_45_DEG_LEFT) && (track_section->flags & (TRACK_DIAGONAL | TRACK_DIAGONAL_2)))end_angle -= 1;
		else if ((track_section->flags & TRACK_EXIT_45_DEG_RIGHT) && !(track_section->flags & (TRACK_DIAGONAL | TRACK_DIAGONAL_2)))end_angle += 1;
		if (end_angle < 0)end_angle += 4;
		if (end_angle > 3)end_angle -= 4;

		int start_tie = start_angle <= 1;
		int end_tie = end_angle > 1;

		//Calculate corrected scale
		double corrected_length = num_meshes * track_type->length;
		if (!start_tie)corrected_length -= track_type->tie_length;
		if (end_tie)corrected_length += track_type->tie_length;
		double corrected_scale = track_section_length / corrected_length;

		double tie_length = corrected_scale * track_type->tie_length;
		double inter_length = corrected_scale * (track_type->length - track_type->tie_length);

		double offset = 0;

		if (extrude_behind)
		{
			num_meshes++;
			offset -= (extrude_behind ? 1 : 0) * corrected_scale * track_type->length;
		}
		if (extrude_in_front)num_meshes++;
		for (int i = 0; i < 2 * num_meshes + 1; i++)
		{
			track_transform_args_t args;
			args.scale = corrected_scale;
			args.offset = offset;
			args.y_offset = 0.0;
			args.z_offset = z_offset;
			args.track_curve = track_section->curve;
			args.flags = track_section->flags;
			args.length = track_section_length;

			const vector2_t uv_offset = vector2(track_mesh.u_offset * i, track_mesh.v_offset * i);

			if ((!(i & 1)) && (i != 0 || start_tie) && (i != 2 * num_meshes || end_tie))
			{
				track_point_t track_point = get_track_point(track_section->curve, track_section->flags, z_offset, args.length, args.offset + track_type->tie_length / 2);
				context_add_model(
					context, &(track_type->tie_mesh),
					transform(
						matrix(track_point.binormal.z, track_point.normal.z, track_point.tangent.z, track_point.binormal.y, track_point.normal.y, track_point.tangent.y, track_point.binormal.x, track_point.normal.x, track_point.tangent.x),
						change_coordinates(track_point.position)
					),
					track_mask
				);
				context_add_model_transformed(context, mesh_tie, track_transform, &args, track_mask, vector2(0.0, 0.0));
				offset += tie_length;
			}
			else if (i & 1)
			{
				int use_alt = i & 2;
				if (track_section->flags & TRACK_ALT_INVERT)use_alt = !use_alt;
				if (extrude_behind)use_alt = !use_alt;
				if (!(track_type->models_loaded & (1 << MODEL_TRACK_ALT)))use_alt = 0;
				//Add track model
				if (!track_mask) {
					if (use_alt)context_add_model_transformed(context, &(track_type->models[MODEL_TRACK_ALT]), track_transform, &args, track_mask, uv_offset);
					else context_add_model_transformed(context, mesh, track_transform, &args, track_mask, uv_offset);
					//Add track mask
				}
				else
				{
					if (start_tie)args.offset = offset - tie_length;
					context_add_model_transformed(context, &(track_type->mask), track_transform, &args, 0, vector2(0.0, 0.0));
				}
				offset += inter_length;
			}
		}




	}
	else
	{
		if (extrude_behind)num_meshes++;
		if (extrude_in_front)num_meshes++;
		for (int i = 0; i < num_meshes; i++)
		{
			track_transform_args_t args;
			args.scale = scale;
			args.offset = (i - (extrude_behind ? 1 : 0)) * length;
			args.y_offset = 0.0;
			args.z_offset = z_offset;
			args.track_curve = track_section->curve;
			args.flags = track_section->flags;
			args.length = track_section_length;

			int alt_available = track_type->models_loaded & (1 << MODEL_TRACK_ALT);
			int use_alt = alt_available && (i & 1);
			if (alt_available && (track_section->flags & TRACK_ALT_INVERT))use_alt = !use_alt;

			const vector2_t uv_offset = vector2(track_mesh.u_offset * i, track_mesh.v_offset * i);

			if (track_mask)
			{
				context_add_model_transformed(context, &(track_type->mask), track_transform, &args, 0, vector2(0.0, 0.0));
			}
			else {
				args.y_offset = track_mesh.y_offset;

				if (use_alt)context_add_model_transformed(context, &(track_type->models[MODEL_TRACK_ALT]), track_transform, &args, track_mask, uv_offset);
				else context_add_model_transformed(context, mesh, track_transform, &args, track_mask, uv_offset);

				if ((track_type->models_loaded & (1 << MODEL_BASE)) && (track_type->flags & TRACK_HAS_SUPPORTS) && !(track_section->flags & TRACK_NO_SUPPORTS))
					context_add_model_transformed(context, &(track_type->models[MODEL_BASE]), base_transform, &args, track_mask, vector2(0.0, 0.0));
				if (track_type->flags & TRACK_SEPARATE_TIE)
				{
					track_point_t track_point = get_track_point(track_section->curve, track_section->flags, z_offset, args.length, args.offset + 0.5 * length);
					context_add_model(
						context, &(track_type->tie_mesh),
						transform_with_distance(
							matrix(track_point.binormal.z, track_point.normal.z, track_point.tangent.z, track_point.binormal.y, track_point.normal.y, track_point.tangent.y, track_point.binormal.x, track_point.normal.x, track_point.tangent.x),
							change_coordinates(track_point.position),
							(args.offset + 0.5 * length) / track_section_length
						),
						track_mask
					);
				}
			}
		}
	}

	if (track_section->flags & TRACK_SPECIAL_MASK)
	{
		int index = get_special_index(track_section->flags);
		if (track_type->models_loaded & (1 << index))
		{
			matrix_t mat = views[1];
			if ((track_section->flags & TRACK_SPECIAL_MASK) != TRACK_SPECIAL_VERTICAL_TWIST_RIGHT && (track_section->flags & TRACK_SPECIAL_MASK) != TRACK_SPECIAL_BARREL_ROLL_RIGHT
				&& (track_section->flags & TRACK_SPECIAL_MASK) != TRACK_SPECIAL_CORKSCREW_RIGHT && (track_section->flags & TRACK_SPECIAL_MASK) != TRACK_SPECIAL_ZERO_G_ROLL_RIGHT
				&& (track_section->flags & TRACK_SPECIAL_MASK) != TRACK_SPECIAL_LARGE_ZERO_G_ROLL_RIGHT)
			{
				mat.entries[6] *= -1;
				mat.entries[7] *= -1;
				mat.entries[8] *= -1;
			}

			if ((track_section->flags & TRACK_SPECIAL_MASK) == TRACK_SPECIAL_BRAKE || (track_section->flags & TRACK_SPECIAL_MASK) == TRACK_SPECIAL_MAGNETIC_BRAKE || (track_section->flags & TRACK_SPECIAL_MASK) == TRACK_SPECIAL_BLOCK_BRAKE || (track_section->flags & TRACK_SPECIAL_MASK) == TRACK_SPECIAL_BOOSTER)
			{
				float special_length = track_type->brake_length;
				if ((track_section->flags & TRACK_SPECIAL_MASK) == TRACK_SPECIAL_BLOCK_BRAKE)special_length = TILE_SIZE;
				int num_special_meshes = (int)floor(0.5 + track_section_length / special_length);
				float special_scale = track_section_length / (num_special_meshes * special_length);
				special_length = special_scale * special_length;
				for (int i = 0; i < num_special_meshes; i++)
				{
					track_transform_args_t args;
					args.scale = special_scale;
					args.offset = i * special_length;
					args.z_offset = z_offset;
					args.track_curve = track_section->curve;
					args.flags = track_section->flags;
					args.length = track_section_length;

					context_add_model_transformed(context, &(track_type->models[index]), track_transform, &args, track_mask, vector2(0.0, 0.0));
				}
			}
			else context_add_model(context, &(track_type->models[index]), transform(mat, vector3(!(track_section->flags & TRACK_VERTICAL) ? -0.5 * TILE_SIZE : 0, z_offset - 2 * CLEARANCE_HEIGHT, 0)), track_mask);
		}
	}

	if ((track_type->flags & TRACK_HAS_SUPPORTS) && !(track_section->flags & TRACK_NO_SUPPORTS))
	{
		int num_supports = (int)floor(0.5 + track_section_length / track_type->support_spacing);
		float support_step = track_section_length / num_supports;
		int entry = 0;
		int exit = 0;
		if (track_section->flags & TRACK_ENTRY_BANK_LEFT)entry = DENOM;
		else if (track_section->flags & TRACK_ENTRY_BANK_RIGHT)entry = -DENOM;
		else entry = 0;

		if (track_section->flags & TRACK_EXIT_BANK_LEFT)exit = DENOM;
		else if (track_section->flags & TRACK_EXIT_BANK_RIGHT)exit = -DENOM;
		else exit = 0;

		for (int i = 0; i < num_supports + 1; i++)
		{
			int u = (i * DENOM) / num_supports;
			int bank_angle = (entry * (DENOM - u) + (exit * u)) / DENOM;

			track_point_t track_point = get_track_point(track_section->curve, track_section->flags, z_offset, track_section_length, i * support_step);

			track_point_t support_point = only_yaw(track_point);

			matrix_t rotation =
				matrix(support_point.binormal.x, support_point.normal.x, support_point.tangent.x, support_point.binormal.y, support_point.normal.y, support_point.tangent.y, support_point.binormal.z, support_point.normal.z, support_point.tangent.z);
			if (bank_angle >= 0)rotation = matrix_mult(views[2], rotation);

			vector3_t translation = change_coordinates(support_point.position);
			translation.y -= track_type->pivot / sqrt(track_point.tangent.x * track_point.tangent.x + track_point.tangent.z * track_point.tangent.z) - track_type->pivot;

			context_add_model(context, &(track_type->models[get_support_index(bank_angle)]), transform(rotation, translation), track_mask);
		}
	}
}

void render_track_section(context_t* context,track_section_t* track_section,track_type_t* track_type,int extrude_behind,int extrude_in_front,int track_mask,int rendered_views,image_t* images,int subtype)
{
	context_begin_render(context);

	for (track_mesh_t& track_mesh : track_type->track_meshes)
	{
		add_track_models(context, track_section, track_type, extrude_behind, extrude_in_front, track_mask, rendered_views, images, subtype, track_mesh);
	}

	context_finalize_render(context);

	for(int i=0; i<4; i++)
	{
		bool fade_shadows = false;
		float fade_shadows_distance = 1.0f;

		for (auto& fade_shadow_entry : track_type->fade_shadow_entries)
		{
			if (fade_shadow_entry.track_section_name == track_section->name && fade_shadow_entry.view == i)
			{
				fade_shadows = true;
				fade_shadows_distance = fade_shadow_entry.distance;
				break;
			}
		}

		if(rendered_views&(1<<i))
		{
			if (track_mask)
			{
				context_render_silhouette(context, rotate_y(0.5 * i * M_PI), images + i, track_type->edge_distance);
			}
			else
			{
				context_render_view(
					context,
					rotate_y(0.5 * i * M_PI),
					images + i,
					track_type->edge_distance,
					track_type->remappable_to_grayscale,
					track_type->remappable_to_grayscale_threshold,
					fade_shadows,
					fade_shadows_distance
				);
			}
		}
	}
	context_end_render(context);
}

int is_in_mask(int x, int y, const mask_t& mask)
{
	x += (mask.image.width / 2);
	if (mask.flipped) {
		x = -x - 1;
		y += 1;
		// not sure why this works
	}
	if (mask.image.pixels[(y + (mask.image.height / 2)) * mask.image.width + x] == mask.color) {
		return 1;
	}
	return 0;
}

int compare_vec(vector3_t vec1,vector3_t vec2,int rot)
{
	return vector3_norm(vector3_sub(vec1,vector3_normalize(matrix_vector(views[rot],vec2))))<0.15;
}

int offset_table_index_with_rot(track_point_t track,int rot)
{

	//Get bank angle
	int banked=fabs(fabs(asin(sqrt(track.normal.x*track.normal.x+track.normal.z*track.normal.z)))-0.25*M_PI)<0.1;
	int right=(banked&&track.binormal.y<0) ? 0x10 : 0;
	//printf("%f %f %f\n",track.tangent.x,track.tangent.y,track.tangent.z);
	//Flat
	if(compare_vec(track.tangent,vector3(0,0,TILE_SIZE),rot))
	{
		//Inverted
		if(track.normal.y<-0.9)return 5;
		//Banked
		else if(banked)return right|3;
		//Unbanked
		else return 0;
	}
	//Gentle
	else if(compare_vec(track.tangent,vector3(0,2*CLEARANCE_HEIGHT,TILE_SIZE),rot))
	{
		if(banked)return right|4;
		else return 1;
	}
	//Steep
	else if(compare_vec(track.tangent,vector3(0,8*CLEARANCE_HEIGHT,TILE_SIZE),rot))return 2;
	//Diagonal flat
	else if(compare_vec(track.tangent,vector3(-TILE_SIZE,0,TILE_SIZE),rot))
	{
		if(banked)return right|7;
		return 6;
	}
	//Diagonal gentle
	else if(compare_vec(track.tangent,vector3(-TILE_SIZE,2*CLEARANCE_HEIGHT,TILE_SIZE),rot)&&!banked)return 8;
	//Diagonal gentle bank
	else if (compare_vec(track.tangent, vector3(-TILE_SIZE, 2 * CLEARANCE_HEIGHT, TILE_SIZE), rot) && banked)return right | 9;
	//Inverted
	else if (track.normal.y < -0.9)return 5;
	// Vertical
	else if (compare_vec(track.normal, vector3(0, 0, -1), rot))
		return 10;
	return 0xFF;
}

int offset_table_index(track_point_t track)
{
	//Check straight
	int index=offset_table_index_with_rot(track,0);
	if(index !=0xFF)return index;

	//Check left
	index=offset_table_index_with_rot(track,1);
	if(index !=0xFF)return 0x60|index;

	//Check right
	index=offset_table_index_with_rot(track,3);
	if(index !=0xFF)return 0x20|index;
	return 0xFF;
}

vector3_t get_offset(float offsets[10][8], int table,int view_angle)
{
	int index=table&0xF;
	int end_angle=table>>5;
	int right=(table&0x10)>>4;
	//printf("index %d end angle %d right %d\n",index,end_angle,right);
	int rotated_view_angle=(view_angle+end_angle+2*right) % 4;
	//printf("view %d rotated view %d\n",view_angle,rotated_view_angle);

	vector3_t offset=vector3(0,0,0);
	if(table ==0xFF)return offset;

	offset.x=0;
	offset.z= offsets[index][2*rotated_view_angle]*TILE_SIZE/32.0;
	offset.y= offsets[index][2*rotated_view_angle+1]*CLEARANCE_HEIGHT/8.0;

	//Check if right banked
	if(right)
	{
		offset.z*=-1;
	}

	//Check if diagonal
	if(index >=6&&index <=8)
	{
		offset.z*=M_SQRT1_2;
		offset.x=offset.z;
	}

	//Correct for end rotation
	if(end_angle !=0)offset=matrix_vector(rotate_y(-0.5*M_PI*end_angle),offset);

	//printf("Offset %d %d %.2f\n",index,rotated_view_angle,offset_tables[index][2*rotated_view_angle+1]);
	return offset;
}

void set_offset(float offsets[10][8], int view_angle,track_section_t* track_section)
{
	int start_table=offset_table_index(track_section->curve(0));
	int end_table=offset_table_index(track_section->curve(track_section->length));

	start_offset=get_offset(offsets, start_table,view_angle);
	end_offset=get_offset(offsets, end_table,view_angle);
}

void render_track_sections(context_t* context,track_section_t* track_section,track_type_t* track_type,int track_mask,int subtype,int views,image_t* sprites)
{
int extrude_behind=track_section->flags&TRACK_EXTRUDE_BEHIND;
int extrude_in_front_even=!(track_section->flags&TRACK_EXIT_45_DEG_LEFT)&&(track_section->flags&TRACK_EXTRUDE_IN_FRONT);
int extrude_in_front_odd=(track_section->flags&TRACK_EXIT_45_DEG_LEFT)&&(track_section->flags&TRACK_EXTRUDE_IN_FRONT);

	if(track_type->flags&TRACK_SPECIAL_OFFSETS)
	{
	set_offset(track_type->offsets, 0,track_section);
		if(views&0x1)render_track_section(context,track_section,track_type,extrude_behind,extrude_in_front_even,track_mask,0x1,sprites,subtype);
	set_offset(track_type->offsets, 1,track_section);
		if(views&0x2)render_track_section(context,track_section,track_type,0,extrude_in_front_odd,track_mask,0x2,sprites,subtype);
	set_offset(track_type->offsets, 2,track_section);
		if(views&0x4)render_track_section(context,track_section,track_type,extrude_behind,extrude_in_front_even,track_mask,0x4,sprites,subtype);
	set_offset(track_type->offsets, 3,track_section);
		if(views&0x8)render_track_section(context,track_section,track_type,0,extrude_in_front_odd,track_mask,0x8,sprites,subtype);
	return;
	}

	if((track_section->flags&TRACK_EXTRUDE_BEHIND)||(track_section->flags&TRACK_EXTRUDE_IN_FRONT))
	{
		if(track_type->flags&TRACK_SEPARATE_TIE)
		{
			if(views&0x1)render_track_section(context,track_section,track_type,extrude_behind,extrude_in_front_even,track_mask,0x1,sprites,subtype);
			if(views&0x2)render_track_section(context,track_section,track_type,0,extrude_in_front_odd,track_mask,0x2,sprites,subtype);
			if(views&0x4)render_track_section(context,track_section,track_type,extrude_behind,extrude_in_front_even,track_mask,0x4,sprites,subtype);
			if(views&0x8)render_track_section(context,track_section,track_type,0,extrude_in_front_odd,track_mask,0x8,sprites,subtype);
		}
		else
		{
			if(views&0x5)render_track_section(context,track_section,track_type,extrude_behind,extrude_in_front_even,track_mask,views&0x5,sprites,subtype);
			if(views&0xA)render_track_section(context,track_section,track_type,0,extrude_in_front_odd,track_mask,views&0xA,sprites,subtype);
		}
	}
	else
	{
		if((track_type->flags&TRACK_SEPARATE_TIE)&&(track_section->flags&TRACK_EXIT_90_DEG))
		{
			if(views&0x1)render_track_section(context,track_section,track_type,0,0,track_mask,0x1,sprites,subtype);
			if(views&0x2)render_track_section(context,track_section,track_type,0,0,track_mask,0x2,sprites,subtype);
			if(views&0x4)render_track_section(context,track_section,track_type,0,0,track_mask,0x4,sprites,subtype);
			if(views&0x8)render_track_section(context,track_section,track_type,0,0,track_mask,0x8,sprites,subtype);
		}
		else if((track_type->flags&TRACK_SEPARATE_TIE))
		{
			if(views&0x3)render_track_section(context,track_section,track_type,0,0,track_mask,views&0x3,sprites,subtype);
			if(views&0xC)render_track_section(context,track_section,track_type,0,0,track_mask,views&0xC,sprites,subtype);
		}
		else
		{
			render_track_section(context,track_section,track_type,0,0,track_mask,views,sprites,subtype);
		}
	}
}

std::array<view_t, 4> load_views(const char* track_section_name, const char* mask_directory, const char* track_mask_directory, const bool flip_view) {
	std::array<view_t, 4> views = { view_t {}, view_t {}, view_t {}, view_t {} };
	
	bool use_track_mask_directory = true;
	json_error_t error;
	json_t* views_description = json_load_file(std::format("{}{}.json", track_mask_directory, track_section_name).c_str(), 0, &error);
	if(views_description == nullptr)
	{
		views_description = json_load_file(std::format("{}{}.json", mask_directory, track_section_name).c_str(), 0, &error);
		if (views_description == nullptr)
		{
			printf("Error: %s at line %d column %d\n", error.text, error.line, error.column);
			std::exit(1);
		}
		use_track_mask_directory = false;
	}

	json_t* alternate_track_section_name = json_object_get(views_description, "alternate_track_section_name");
	if (alternate_track_section_name == nullptr || !json_is_string(alternate_track_section_name))
	{
		printf("Error: Property \"alternate_track_section_name\" not found or is not a string\n");
		std::abort();
	}
	const char* alternate_track_section_name_value = json_string_value(alternate_track_section_name);

	json_t* alternate_track_section_flipped = json_object_get(views_description, "alternate_track_section_flipped");
	if (alternate_track_section_flipped == nullptr || !json_is_boolean(alternate_track_section_flipped))
	{
		printf("Error: Property \"alternate_track_section_flipped\" not found or is not a boolean\n");
		std::abort();
	}
	const bool alternate_track_section_flipped_value = json_boolean_value(alternate_track_section_flipped);

	if (strlen(alternate_track_section_name_value) != 0) {
		return load_views(alternate_track_section_name_value, mask_directory, track_mask_directory, alternate_track_section_flipped_value);
	}

	json_t* views_array = json_object_get(views_description, "views");
	if (views_array == nullptr || !json_is_array(views_array))
	{
		printf("Error: Property \"views\" not found or is not an array\n");
		std::abort();
	}
	const int views_array_size = json_array_size(views_array);
	if (views_array_size != 4)
	{
		printf("Error: Property \"views\" array size is not 4\n");
		std::abort();
	}

	for (int i = 0; i < 4; i++) {
		view_t& view = views.at(i);

		json_t* view_description = flip_view ? json_array_get(views_array, 3 - i) : json_array_get(views_array, i);

		json_t* uses_track_mask = json_object_get(view_description, "uses_track_mask");
		if (uses_track_mask == nullptr || !json_is_boolean(uses_track_mask))
		{
			printf("Error: Property \"uses_track_mask\" not found or is not a boolean\n");
			std::abort();
		}
		const bool uses_track_mask_value = json_boolean_value(uses_track_mask);

		json_t* enforce_non_overlapping = json_object_get(view_description, "enforce_non_overlapping");
		if (enforce_non_overlapping == nullptr || !json_is_boolean(enforce_non_overlapping))
		{
			printf("Error: Property \"enforce_non_overlapping\" not found or is not a boolean\n");
			std::abort();
		}
		const bool enforce_non_overlapping_value = json_boolean_value(enforce_non_overlapping);

		if (uses_track_mask_value) {
			view.flags |= VIEW_NEEDS_TRACK_MASK;
		}

		if (enforce_non_overlapping_value) {
			view.flags |= VIEW_ENFORCE_NON_OVERLAPPING;
		}

		json_t* mask_array = json_object_get(view_description, "masks");
		if (mask_array == nullptr || !json_is_array(mask_array))
		{
			printf("Error: Property \"mask_array\" not found or is not an array\n");
			std::abort();
		}
		const int masks_array_size = json_array_size(mask_array);

		for (int mask_i = 0; mask_i < masks_array_size; mask_i++) {
			json_t* mask_description = json_array_get(mask_array, mask_i);

			json_t* use_mask = json_object_get(mask_description, "use_mask");
			if (use_mask == nullptr || !json_is_boolean(use_mask))
			{
				printf("Error: Property \"use_mask\" not found or is not a boolean\n");
				std::abort();
			}
			const bool use_mask_value = json_boolean_value(use_mask);

			if (use_mask_value) {
				json_t* image_path = json_object_get(mask_description, "image_path");
				if (image_path == nullptr || !json_is_string(image_path))
				{
					printf("Error: Property \"image_path\" not found or is not a string\n");
					std::abort();
				}
				const char* image_path_value = json_string_value(image_path);

				json_t* color = json_object_get(mask_description, "color");
				if (color == nullptr || !json_is_integer(color))
				{
					printf("Error: Property \"color\" not found or is not an integer\n");
					std::abort();
				}
				const int color_value = json_integer_value(color);

				json_t* mask_type = json_object_get(mask_description, "mask_type");
				if (mask_type == nullptr || !json_is_string(mask_type))
				{
					printf("Error: Property \"mask_type\" not found or is not a string\n");
					std::abort();
				}
				const char* mask_type_value = json_string_value(mask_type);

				int track_mask_op = TRACK_MASK_NONE;
				if (!strcmp(mask_type_value, "difference")) {
					track_mask_op = TRACK_MASK_DIFFERENCE;
				}
				if (!strcmp(mask_type_value, "intersect")) {
					track_mask_op = TRACK_MASK_INTERSECT;
				}
				if (!strcmp(mask_type_value, "union")) {
					track_mask_op = TRACK_MASK_UNION;
				}
				if (!strcmp(mask_type_value, "transfer_next")) {
					track_mask_op = TRACK_MASK_TRANSFER_NEXT;
				}

				json_t* offset_x = json_object_get(mask_description, "offset_x");
				if (offset_x == nullptr || !json_is_integer(offset_x))
				{
					printf("Error: Property \"offset_x\" not found or is not an integer\n");
					std::abort();
				}
				const int offset_x_value = json_integer_value(offset_x);

				json_t* offset_y = json_object_get(mask_description, "offset_y");
				if (offset_y == nullptr || !json_is_integer(offset_y))
				{
					printf("Error: Property \"offset_y\" not found or is not an integer\n");
					std::abort();
				}
				const int offset_y_value = json_integer_value(offset_y);

				json_t* flipped = json_object_get(mask_description, "flipped");
				if (flipped == nullptr || !json_is_boolean(flipped))
				{
					printf("Error: Property \"flipped\" not found or is not a boolean\n");
					std::abort();
				}
				const bool flipped_value = json_boolean_value(flipped);

				auto image_file_path = std::string();
				if (use_track_mask_directory) {
					image_file_path = std::format("{}{}", track_mask_directory, image_path_value);
				}
				else {
					image_file_path = std::format("{}{}", mask_directory, image_path_value);
				}
				
				FILE* image_file;
				image_file = fopen(image_file_path.c_str(), "rb");
				if (image_file == nullptr) {
					printf("Error: Could not open %s\n", image_file_path.c_str());
					std::abort();
				}

				image_t image;
				if (image_read_png(&image, image_file) != 0) {
					printf("Error: Could not read %s as png\n", image_file_path.c_str());
					std::abort();
				}

				fclose(image_file);

				view.masks.push_back(mask_t{ true, image, uint8_t(color_value), track_mask_op, flip_view ? -offset_x_value : offset_x_value, offset_y_value, flip_view || flipped_value });
			}
			else {
				view.masks.push_back(mask_t{ false, image_t{}, 0, 0, 0, 0, false });
			}
		}
	}

	return views;
}

void write_track_section(context_t* context,track_section_t* track_section,track_type_t* track_type,const char* base_directory,const char* mask_directory,const char* filename,json_t* sprites,int subtype,image_t* overlay)
{
	int z_offset=(int)(track_type->z_offset+0.499999);
	image_t full_sprites[4];
	render_track_sections(context,track_section,track_type,0,subtype,0xF,full_sprites);

	if(overlay !=NULL&&!(track_type->flags&TRACK_NO_LIFT_SPRITE))
	{
		for(int i=0; i<4; i++)image_blit(full_sprites+i,overlay+i,0,track_type->lift_offset-z_offset);
	}

	const std::string track_mask_directory = std::format("{}masks/", base_directory);

	const auto views = load_views(track_section->name, mask_directory, track_mask_directory.c_str(), false);

	image_t track_masks[4];
	int track_mask_views=0;
	for(int i=0; i<4; i++)track_mask_views|=(views[i].flags&VIEW_NEEDS_TRACK_MASK ? 1 : 0)<<i;
	if(track_mask_views !=0)render_track_sections(context,track_section,track_type,1,subtype,track_mask_views,track_masks);

	for(int angle=0; angle<4; angle++)
	{
		if(views[angle].masks.size() == 0) continue;

		const view_t& view = views[angle];

		char final_filename[512];
		char relative_filename[512];
		snprintf(relative_filename,512,"%s_%d.png",filename,angle+1);

		for(int sprite=0; sprite<view.masks.size(); sprite++)
		{
			char final_filename[512];
			char relative_filename[512];
			if(view.masks.size() ==1)snprintf(relative_filename,512,"%s_%d.png",filename,angle+1);
			else snprintf(relative_filename,512,"%s_%d_%d.png",filename,angle+1,sprite+1);
			snprintf(final_filename,512,"%s%s",base_directory,relative_filename);
			//y		snprintf(final_filename,512,"../ImageEncode/%s",relative_filename);
			printf("%s\n",final_filename);

			image_t part_sprite;
			image_copy(full_sprites+angle,&part_sprite);

			const mask_t& mask = view.masks.at(sprite);

			if(mask.use_mask)
			{

				for(int x=0; x<full_sprites[angle].width; x++)
					for(int y=0; y<full_sprites[angle].height; y++)
					{
						int in_mask=is_in_mask(x+full_sprites[angle].x_offset,y+full_sprites[angle].y_offset+((track_section->flags&TRACK_OFFSET_SPRITE_MASK) ? (z_offset-8) : 0), mask);

						if(mask.track_mask_op !=TRACK_MASK_NONE)
						{
							int mask_x=(x+full_sprites[angle].x_offset)-track_masks[angle].x_offset;
							int mask_y=(y+full_sprites[angle].y_offset)-track_masks[angle].y_offset;

							int in_track_mask = 0;
							if (mask_x >= 0 && mask_y >= 0 && mask_x < track_masks[angle].width && mask_y < track_masks[angle].height) {
								const float track_depth = full_sprites[angle].depths[x + (y * full_sprites[angle].width)];
								const float mask_depth = track_masks[angle].depths[mask_x + (mask_y * track_masks[angle].width)];

								in_track_mask = mask_depth > 0.0 && mask_depth <= track_depth ? 1 : 0;
							}
							
							switch (mask.track_mask_op)
							{
							case TRACK_MASK_DIFFERENCE:
								in_mask=in_mask&&!in_track_mask;
								break;
							case TRACK_MASK_INTERSECT:
								in_mask=in_mask&&in_track_mask;
								break;
							case TRACK_MASK_UNION:
								in_mask=in_mask||in_track_mask;
								break;
							}

							if(sprite < (view.masks.size() - 1) && (mask.track_mask_op & TRACK_MASK_TRANSFER_NEXT) && in_track_mask
							  &&is_in_mask(x+full_sprites[angle].x_offset,y+full_sprites[angle].y_offset+((track_section->flags&TRACK_OFFSET_SPRITE_MASK) ? (z_offset-8) : 0), view.masks.at(sprite+1)))
								in_mask=1;
						}

						if(view.flags&VIEW_ENFORCE_NON_OVERLAPPING)
						{
							for(int i=0; i<sprite; i++)
							{
								if(is_in_mask(x+full_sprites[angle].x_offset,y+full_sprites[angle].y_offset+((track_section->flags&TRACK_OFFSET_SPRITE_MASK) ? (z_offset-8) : 0), view.masks.at(i)))in_mask=0;//Note z offset untested
							}
						}

						if(in_mask)
						{
							part_sprite.pixels[x+part_sprite.width*y]=full_sprites[angle].pixels[x+full_sprites[angle].width*y];
						}
						else
						{
							part_sprite.pixels[x+part_sprite.width*y]=0;
						}
					}
				part_sprite.x_offset += mask.x_offset;
				part_sprite.y_offset += mask.y_offset;
			}

			FILE* file=fopen(final_filename,"wb");
			if(file ==NULL)
			{
				printf("Error: could not open %s for writing\n",final_filename);
				exit(1);
			}
			//if(view->flags&VIEW_NEEDS_TRACK_MASK)image_write_png(&(track_masks[angle]),file);
			image_crop(&part_sprite);
			image_write_png(&part_sprite,file);
			//image_write_png(full_sprites+angle,file);
			fclose(file);

			json_t* sprite_entry=json_object();
			json_object_set(sprite_entry,"path",json_string(relative_filename));
			json_object_set(sprite_entry,"x",json_integer(part_sprite.x_offset));
			json_object_set(sprite_entry,"y",json_integer(part_sprite.y_offset));
			json_object_set(sprite_entry,"palette",json_string("keep"));
			json_array_append(sprites,sprite_entry);
			image_destroy(&part_sprite);
		}

		if(view.flags&VIEW_NEEDS_TRACK_MASK)image_destroy(track_masks+angle);
		image_destroy(full_sprites+angle);
	}
}

int write_track_subtype(context_t* context,track_type_t* track_type,track_list_t track_list,json_t* sprites,const char* base_dir,const char* mask_dir,const char* output_dir,int subtype)
{
char output_path[300];
const char* suffix="";

uint64_t groups=0;
	switch(subtype)
	{
	case TRACK_SUBTYPE_DEFAULT:
		groups=track_type->groups;
		suffix="";
		break;
	case TRACK_SUBTYPE_LIFT:
		groups=track_type->lift_groups;
		suffix="_lift";
		break;
	}

	//Flat
	if(groups&TRACK_GROUP_FLAT)
	{
	sprintf(output_path,"%.255sflat%s",output_dir,suffix);
	if(subtype ==TRACK_SUBTYPE_LIFT)write_track_section(context,&(track_list.flat_asymmetric),track_type,base_dir,mask_dir,output_path,sprites,subtype,flat_chain);
	else write_track_section(context,&(track_list.flat),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	}
	if(groups&TRACK_GROUP_BRAKES)
	{
	sprintf(output_path,"%.255sbrake%s",output_dir,suffix);
	write_track_section(context,&(track_list.brake),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	}
	if(groups&TRACK_GROUP_BLOCK_BRAKES)
	{
	sprintf(output_path,"%.255sblock_brake%s",output_dir,suffix);
	write_track_section(context,&(track_list.block_brake),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	}
	if(groups&TRACK_GROUP_SLOPED_BRAKES)
	{
	sprintf(output_path,"%.255sbrake_gentle%s",output_dir,suffix);
	write_track_section(context,&(track_list.brake_gentle),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	}
	if(groups&TRACK_GROUP_MAGNETIC_BRAKES)
	{
	sprintf(output_path,"%.255smagnetic_brake%s",output_dir,suffix);
	write_track_section(context,&(track_list.magnetic_brake),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	}


	if(groups&TRACK_GROUP_BOOSTERS)
	{
	sprintf(output_path,"%.255sbooster%s",output_dir,suffix);
	write_track_section(context,&(track_list.booster),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	}
	//Launched lift
	if(groups&TRACK_GROUP_LAUNCHED_LIFTS)
	{
	sprintf(output_path,"%.255spowered_lift%s",output_dir,suffix);
	write_track_section(context,&(track_list.launched_lift),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	}
	if(groups&TRACK_GROUP_VERTICAL_BOOSTERS)
	{
	sprintf(output_path,"%.255svertical_booster%s",output_dir,suffix);
	write_track_section(context,&(track_list.vertical_booster),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	}

	//Slopes
	if(groups&TRACK_GROUP_GENTLE_SLOPES)
	{
	sprintf(output_path,"%.255sflat_to_gentle_up%s",output_dir,suffix);
	write_track_section(context,&(track_list.flat_to_gentle_up),track_type,base_dir,mask_dir,output_path,sprites,subtype,subtype ==TRACK_SUBTYPE_LIFT ? flat_to_gentle_up_chain : NULL);
	sprintf(output_path,"%.255sgentle_up_to_flat%s",output_dir,suffix);
	write_track_section(context,&(track_list.gentle_up_to_flat),track_type,base_dir,mask_dir,output_path,sprites,subtype,subtype ==TRACK_SUBTYPE_LIFT ? gentle_up_to_flat_chain : NULL);
	sprintf(output_path,"%.255sgentle%s",output_dir,suffix);
	write_track_section(context,&(track_list.gentle),track_type,base_dir,mask_dir,output_path,sprites,subtype,subtype ==TRACK_SUBTYPE_LIFT ? gentle_chain : NULL);
	}
	//TODO should probably be inside slopes
		if(groups&TRACK_GROUP_MAGNETIC_BRAKES)
		{
		sprintf(output_path,"%.255smagnetic_brake_gentle%s",output_dir,suffix);
		write_track_section(context,&(track_list.magnetic_brake_gentle),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
		}

	if(groups&TRACK_GROUP_STEEP_SLOPES)
	{
	sprintf(output_path,"%.255sgentle_to_steep_up%s",output_dir,suffix);
	write_track_section(context,&(track_list.gentle_to_steep_up),track_type,base_dir,mask_dir,output_path,sprites,subtype,subtype ==TRACK_SUBTYPE_LIFT ? gentle_to_steep_up_chain : NULL);
	sprintf(output_path,"%.255ssteep_to_gentle_up%s",output_dir,suffix);
	write_track_section(context,&(track_list.steep_to_gentle_up),track_type,base_dir,mask_dir,output_path,sprites,subtype,subtype ==TRACK_SUBTYPE_LIFT ? steep_to_gentle_up_chain : NULL);
	sprintf(output_path,"%.255ssteep%s",output_dir,suffix);
	write_track_section(context,&(track_list.steep),track_type,base_dir,mask_dir,output_path,sprites,subtype,subtype ==TRACK_SUBTYPE_LIFT ? steep_chain : NULL);
	}
	
	if(groups&TRACK_GROUP_VERTICAL_SLOPES)
	{
	sprintf(output_path,"%.255ssteep_to_vertical_up%s",output_dir,suffix);
	write_track_section(context,&(track_list.steep_to_vertical_up),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	sprintf(output_path,"%.255svertical_to_steep_up%s",output_dir,suffix);
	write_track_section(context,&(track_list.vertical_to_steep_up),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	sprintf(output_path,"%.255svertical%s",output_dir,suffix);
	write_track_section(context,&(track_list.vertical),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	}

	//Turns
	if(groups&TRACK_GROUP_TURNS)
	{
	sprintf(output_path,"%.255ssmall_turn_left%s",output_dir,suffix);
	write_track_section(context,&(track_list.small_turn_left),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	sprintf(output_path,"%.255smedium_turn_left%s",output_dir,suffix);
	write_track_section(context,&(track_list.medium_turn_left),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	sprintf(output_path,"%.255slarge_turn_left_to_diag%s",output_dir,suffix);
	write_track_section(context,&(track_list.large_turn_left_to_diag),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	sprintf(output_path,"%.255slarge_turn_right_to_diag%s",output_dir,suffix);
	write_track_section(context,&(track_list.large_turn_right_to_diag),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	}

	//Diagonals
	if(groups&TRACK_GROUP_DIAGONALS)
	{
	sprintf(output_path,"%.255sflat_diag%s",output_dir,suffix);
	write_track_section(context,&(track_list.flat_diag),track_type,base_dir,mask_dir,output_path,sprites,subtype,subtype ==TRACK_SUBTYPE_LIFT ? flat_diag_chain : NULL);
	}
	if(groups&TRACK_GROUP_DIAGONAL_BRAKES)
	{
		if(groups&TRACK_GROUP_BRAKES)
		{
		sprintf(output_path,"%.255sbrake_diag%s",output_dir,suffix);
		write_track_section(context,&(track_list.brake_diag),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
		}
		if(groups&TRACK_GROUP_BLOCK_BRAKES)
		{
		sprintf(output_path,"%.255sblock_brake_diag%s",output_dir,suffix);
		write_track_section(context,&(track_list.block_brake_diag),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
		}
		if(groups&TRACK_GROUP_MAGNETIC_BRAKES)
		{
		sprintf(output_path,"%.255smagnetic_brake_diag%s",output_dir,suffix);
		write_track_section(context,&(track_list.magnetic_brake_diag),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
		}
	};
	if((groups&TRACK_GROUP_DIAGONALS)&&(groups&TRACK_GROUP_GENTLE_SLOPES))
	{
	sprintf(output_path,"%.255sflat_to_gentle_up_diag%s",output_dir,suffix);
	write_track_section(context,&(track_list.flat_to_gentle_up_diag),track_type,base_dir,mask_dir,output_path,sprites,subtype,subtype ==TRACK_SUBTYPE_LIFT ? flat_to_gentle_up_diag_chain : NULL);
	sprintf(output_path,"%.255sgentle_to_flat_up_diag%s",output_dir,suffix);
	write_track_section(context,&(track_list.gentle_to_flat_up_diag),track_type,base_dir,mask_dir,output_path,sprites,subtype,subtype ==TRACK_SUBTYPE_LIFT ? gentle_to_flat_up_diag_chain : NULL);
	sprintf(output_path,"%.255sgentle_diag%s",output_dir,suffix);
	write_track_section(context,&(track_list.gentle_diag),track_type,base_dir,mask_dir,output_path,sprites,subtype,subtype ==TRACK_SUBTYPE_LIFT ? gentle_diag_chain : NULL);
	};
	if(groups&TRACK_GROUP_DIAGONAL_BRAKES)
	{
		if(groups&TRACK_GROUP_SLOPED_BRAKES)
		{
		sprintf(output_path,"%.255sbrake_gentle_diag%s",output_dir,suffix);
		write_track_section(context,&(track_list.brake_gentle_diag),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
		}
		if(groups&TRACK_GROUP_MAGNETIC_BRAKES)
		{
		sprintf(output_path,"%.255smagnetic_brake_gentle_diag%s",output_dir,suffix);
		write_track_section(context,&(track_list.magnetic_brake_gentle_diag),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
		}
	};
	if((groups&TRACK_GROUP_DIAGONALS)&&(groups&TRACK_GROUP_STEEP_SLOPES))
	{
	sprintf(output_path,"%.255sgentle_to_steep_up_diag%s",output_dir,suffix);
	write_track_section(context,&(track_list.gentle_to_steep_up_diag),track_type,base_dir,mask_dir,output_path,sprites,subtype,subtype ==TRACK_SUBTYPE_LIFT ? gentle_to_steep_up_diag_chain : NULL);
	sprintf(output_path,"%.255ssteep_to_gentle_up_diag%s",output_dir,suffix);
	write_track_section(context,&(track_list.steep_to_gentle_up_diag),track_type,base_dir,mask_dir,output_path,sprites,subtype,subtype ==TRACK_SUBTYPE_LIFT ? steep_to_gentle_up_diag_chain : NULL);
	sprintf(output_path,"%.255ssteep_diag%s",output_dir,suffix);
	write_track_section(context,&(track_list.steep_diag),track_type,base_dir,mask_dir,output_path,sprites,subtype,subtype ==TRACK_SUBTYPE_LIFT ? steep_diag_chain : NULL);
	}


/*
	if((groups&TRACK_GROUP_DIAGONALS)&&(groups&TRACK_GROUP_VERTICAL_SLOPES))
	{
	sprintf(output_path,"%.255ssteep_to_vertical_up_diag%s",output_dir,suffix);
	write_track_section(context,&(track_list.steep_to_vertical_up_diag),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	sprintf(output_path,"%.255svertical_to_steep_up_diag%s",output_dir,suffix);
	write_track_section(context,&(track_list.vertical_to_steep_up_diag),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	sprintf(output_path,"%.255svertical_diag%s",output_dir,suffix);
	write_track_section(context,&(track_list.vertical),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	sprintf(output_path,"%.255svertical_twist_left_to_diag_up%s",output_dir,suffix);
	write_track_section(context,&(track_list.vertical_twist_left_to_diag_up),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	sprintf(output_path,"%.255svertical_twist_right_to_diag_up%s",output_dir,suffix);
	write_track_section(context,&(track_list.vertical_twist_right_to_diag_up),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	sprintf(output_path,"%.255svertical_twist_left_to_orthogonal_up%s",output_dir,suffix);
	write_track_section(context,&(track_list.vertical_twist_left_to_orthogonal_up),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	sprintf(output_path,"%.255svertical_twist_right_to_orthogonal_up%s",output_dir,suffix);
	write_track_section(context,&(track_list.vertical_twist_right_to_orthogonal_up),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	}
*/
	//Banked turns
	if(groups&TRACK_GROUP_BANKED_TURNS)
	{
	sprintf(output_path,"%.255sflat_to_left_bank%s",output_dir,suffix);
	write_track_section(context,&(track_list.flat_to_left_bank),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	sprintf(output_path,"%.255sflat_to_right_bank%s",output_dir,suffix);
	write_track_section(context,&(track_list.flat_to_right_bank),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	sprintf(output_path,"%.255sleft_bank_to_gentle_up%s",output_dir,suffix);
	write_track_section(context,&(track_list.left_bank_to_gentle_up),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	sprintf(output_path,"%.255sright_bank_to_gentle_up%s",output_dir,suffix);
	write_track_section(context,&(track_list.right_bank_to_gentle_up),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	sprintf(output_path,"%.255sgentle_up_to_left_bank%s",output_dir,suffix);
	write_track_section(context,&(track_list.gentle_up_to_left_bank),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	sprintf(output_path,"%.255sgentle_up_to_right_bank%s",output_dir,suffix);
	write_track_section(context,&(track_list.gentle_up_to_right_bank),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);

	sprintf(output_path,"%.255sleft_bank%s",output_dir,suffix);
	write_track_section(context,&(track_list.left_bank),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);

		if(groups&TRACK_GROUP_DIAGONALS)
		{
		sprintf(output_path,"%.255sflat_to_left_bank_diag%s",output_dir,suffix);
		write_track_section(context,&(track_list.flat_to_left_bank_diag),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
		sprintf(output_path,"%.255sflat_to_right_bank_diag%s",output_dir,suffix);
		write_track_section(context,&(track_list.flat_to_right_bank_diag),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
		sprintf(output_path,"%.255sleft_bank_to_gentle_up_diag%s",output_dir,suffix);
		write_track_section(context,&(track_list.left_bank_to_gentle_up_diag),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
		sprintf(output_path,"%.255sright_bank_to_gentle_up_diag%s",output_dir,suffix);
		write_track_section(context,&(track_list.right_bank_to_gentle_up_diag),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
		sprintf(output_path,"%.255sgentle_up_to_left_bank_diag%s",output_dir,suffix);
		write_track_section(context,&(track_list.gentle_up_to_left_bank_diag),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
		sprintf(output_path,"%.255sgentle_up_to_right_bank_diag%s",output_dir,suffix);
		write_track_section(context,&(track_list.gentle_up_to_right_bank_diag),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
		sprintf(output_path,"%.255sleft_bank_diag%s",output_dir,suffix);
		write_track_section(context,&(track_list.left_bank_diag),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
		}

	sprintf(output_path,"%.255ssmall_turn_left_bank%s",output_dir,suffix);
	write_track_section(context,&(track_list.small_turn_left_bank),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	sprintf(output_path,"%.255smedium_turn_left_bank%s",output_dir,suffix);
	write_track_section(context,&(track_list.medium_turn_left_bank),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	sprintf(output_path,"%.255slarge_turn_left_to_diag_bank%s",output_dir,suffix);
	write_track_section(context,&(track_list.large_turn_left_to_diag_bank),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	sprintf(output_path,"%.255slarge_turn_right_to_diag_bank%s",output_dir,suffix);
	write_track_section(context,&(track_list.large_turn_right_to_diag_bank),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	}

	//Sloped turns
	if(groups&TRACK_GROUP_SLOPED_TURNS&&(groups&TRACK_GROUP_GENTLE_SLOPES))
	{
	sprintf(output_path,"%.255ssmall_turn_left_gentle_up%s",output_dir,suffix);
	write_track_section(context,&(track_list.small_turn_left_gentle_up),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	sprintf(output_path,"%.255ssmall_turn_right_gentle_up%s",output_dir,suffix);
	write_track_section(context,&(track_list.small_turn_right_gentle_up),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	sprintf(output_path,"%.255smedium_turn_left_gentle_up%s",output_dir,suffix);
	write_track_section(context,&(track_list.medium_turn_left_gentle_up),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	sprintf(output_path,"%.255smedium_turn_right_gentle_up%s",output_dir,suffix);
	write_track_section(context,&(track_list.medium_turn_right_gentle_up),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	}
	if((groups&TRACK_GROUP_STEEP_SLOPED_TURNS)&&(groups&TRACK_GROUP_STEEP_SLOPES))
	{
	sprintf(output_path,"%.255svery_small_turn_left_steep_up%s",output_dir,suffix);
	write_track_section(context,&(track_list.very_small_turn_left_steep_up),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	sprintf(output_path,"%.255svery_small_turn_right_steep_up%s",output_dir,suffix);
	write_track_section(context,&(track_list.very_small_turn_right_steep_up),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	}
	if((groups&TRACK_GROUP_SLOPED_TURNS)&&(groups&TRACK_GROUP_VERTICAL_SLOPES))
	{
	sprintf(output_path,"%.255svertical_twist_left_up%s",output_dir,suffix);
	write_track_section(context,&(track_list.vertical_twist_left_up),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	sprintf(output_path,"%.255svertical_twist_right_up%s",output_dir,suffix);
	write_track_section(context,&(track_list.vertical_twist_right_up),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	}

	//Sloped banked turns

	if(groups&TRACK_GROUP_BANKED_SLOPED_TURNS)
	{
	sprintf(output_path,"%.255sgentle_up_to_gentle_up_left_bank%s",output_dir,suffix);
	write_track_section(context,&(track_list.gentle_up_to_gentle_up_left_bank),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	sprintf(output_path,"%.255sgentle_up_to_gentle_up_right_bank%s",output_dir,suffix);
	write_track_section(context,&(track_list.gentle_up_to_gentle_up_right_bank),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	sprintf(output_path,"%.255sgentle_up_left_bank_to_gentle_up%s",output_dir,suffix);
	write_track_section(context,&(track_list.gentle_up_left_bank_to_gentle_up),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	sprintf(output_path,"%.255sgentle_up_right_bank_to_gentle_up%s",output_dir,suffix);
	write_track_section(context,&(track_list.gentle_up_right_bank_to_gentle_up),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	sprintf(output_path,"%.255sleft_bank_to_gentle_up_left_bank%s",output_dir,suffix);
	write_track_section(context,&(track_list.left_bank_to_gentle_up_left_bank),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	sprintf(output_path,"%.255sright_bank_to_gentle_up_right_bank%s",output_dir,suffix);
	write_track_section(context,&(track_list.right_bank_to_gentle_up_right_bank),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	sprintf(output_path,"%.255sgentle_up_left_bank_to_left_bank%s",output_dir,suffix);
	write_track_section(context,&(track_list.gentle_up_left_bank_to_left_bank),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	sprintf(output_path,"%.255sgentle_up_right_bank_to_right_bank%s",output_dir,suffix);
	write_track_section(context,&(track_list.gentle_up_right_bank_to_right_bank),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	sprintf(output_path,"%.255sgentle_up_left_bank%s",output_dir,suffix);
	write_track_section(context,&(track_list.gentle_up_left_bank),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	sprintf(output_path,"%.255sgentle_up_right_bank%s",output_dir,suffix);
	write_track_section(context,&(track_list.gentle_up_right_bank),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	sprintf(output_path,"%.255sflat_to_gentle_up_left_bank%s",output_dir,suffix);
	write_track_section(context,&(track_list.flat_to_gentle_up_left_bank),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	sprintf(output_path,"%.255sflat_to_gentle_up_right_bank%s",output_dir,suffix);
	write_track_section(context,&(track_list.flat_to_gentle_up_right_bank),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	sprintf(output_path,"%.255sgentle_up_left_bank_to_flat%s",output_dir,suffix);
	write_track_section(context,&(track_list.gentle_up_left_bank_to_flat),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	sprintf(output_path,"%.255sgentle_up_right_bank_to_flat%s",output_dir,suffix);
	write_track_section(context,&(track_list.gentle_up_right_bank_to_flat),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);

	sprintf(output_path,"%.255ssmall_turn_left_bank_gentle_up%s",output_dir,suffix);
	write_track_section(context,&(track_list.small_turn_left_bank_gentle_up),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	sprintf(output_path,"%.255ssmall_turn_right_bank_gentle_up%s",output_dir,suffix);
	write_track_section(context,&(track_list.small_turn_right_bank_gentle_up),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	sprintf(output_path,"%.255smedium_turn_left_bank_gentle_up%s",output_dir,suffix);
	write_track_section(context,&(track_list.medium_turn_left_bank_gentle_up),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	sprintf(output_path,"%.255smedium_turn_right_bank_gentle_up%s",output_dir,suffix);
	write_track_section(context,&(track_list.medium_turn_right_bank_gentle_up),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	}

	//Miscellaneous
	if(groups&TRACK_GROUP_S_BENDS)
	{
		sprintf(output_path,"%.255ss_bend_left%s",output_dir,suffix);
		write_track_section(context,&(track_list.s_bend_left),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
		sprintf(output_path,"%.255ss_bend_right%s",output_dir,suffix);
		write_track_section(context,&(track_list.s_bend_right),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	}
	if(groups&TRACK_GROUP_BANKED_S_BENDS)
	{
		sprintf(output_path,"%.255ss_bend_bank_left%s",output_dir,suffix);
		write_track_section(context,&(track_list.s_bend_left_bank),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
		sprintf(output_path,"%.255ss_bend_bank_right%s",output_dir,suffix);
		write_track_section(context,&(track_list.s_bend_right_bank),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	}

	if(groups&TRACK_GROUP_HELICES)
	{
		sprintf(output_path,"%.255ssmall_helix_left_up%s",output_dir,suffix);
		write_track_section(context,&(track_list.small_helix_left_up),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
		sprintf(output_path,"%.255ssmall_helix_right_up%s",output_dir,suffix);
		write_track_section(context,&(track_list.small_helix_right_up),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
		sprintf(output_path,"%.255smedium_helix_left_up%s",output_dir,suffix);
		write_track_section(context,&(track_list.medium_helix_left_up),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
		sprintf(output_path,"%.255smedium_helix_right_up%s",output_dir,suffix);
		write_track_section(context,&(track_list.medium_helix_right_up),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	}

	//Inversions
	if(groups&TRACK_GROUP_BARREL_ROLLS)
	{
		sprintf(output_path,"%.255sbarrel_roll_left%s",output_dir,suffix);
		write_track_section(context,&(track_list.barrel_roll_left),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
		sprintf(output_path,"%.255sbarrel_roll_right%s",output_dir,suffix);
		write_track_section(context,&(track_list.barrel_roll_right),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	}
	if(groups&TRACK_GROUP_INLINE_TWISTS)
	{
		sprintf(output_path,"%.255sinline_twist_left%s",output_dir,suffix);
		write_track_section(context,&(track_list.inline_twist_left),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
		sprintf(output_path,"%.255sinline_twist_right%s",output_dir,suffix);
		write_track_section(context,&(track_list.inline_twist_right),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	}
	if(groups&TRACK_GROUP_HALF_LOOPS)
	{
		sprintf(output_path,"%.255shalf_loop%s",output_dir,suffix);
		write_track_section(context,&(track_list.half_loop),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	}
	if(groups&TRACK_GROUP_VERTICAL_LOOPS)
	{
		sprintf(output_path,"%.255sleft_vertical_loop%s",output_dir,suffix);
		write_track_section(context,&(track_list.left_vertical_loop),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
		sprintf(output_path,"%.255sright_vertical_loop%s",output_dir,suffix);
		write_track_section(context,&(track_list.right_vertical_loop),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	}
	if(groups&TRACK_GROUP_LARGE_SLOPE_TRANSITIONS)
	{
		sprintf(output_path,"%.255sflat_to_steep_up%s",output_dir,suffix);
		write_track_section(context,&(track_list.flat_to_steep_up),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
		sprintf(output_path,"%.255ssteep_to_flat_up%s",output_dir,suffix);
		write_track_section(context,&(track_list.steep_to_flat_up),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	}
	if(groups&TRACK_GROUP_QUARTER_LOOPS)
	{
		sprintf(output_path,"%.255squarter_loop_up%s",output_dir,suffix);
		write_track_section(context,&(track_list.quarter_loop_up),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	}
	if(groups&TRACK_GROUP_CORKSCREWS)
	{
		sprintf(output_path,"%.255scorkscrew_left%s",output_dir,suffix);
		write_track_section(context,&(track_list.corkscrew_left),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
		sprintf(output_path,"%.255scorkscrew_right%s",output_dir,suffix);
		write_track_section(context,&(track_list.corkscrew_right),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	}
	if(groups&TRACK_GROUP_LARGE_CORKSCREWS)
	{
		sprintf(output_path,"%.255slarge_corkscrew_left%s",output_dir,suffix);
		write_track_section(context,&(track_list.large_corkscrew_left),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
		sprintf(output_path,"%.255slarge_corkscrew_right%s",output_dir,suffix);
		write_track_section(context,&(track_list.large_corkscrew_right),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	}
	if(groups&TRACK_GROUP_TURN_BANK_TRANSITIONS)
	{
		sprintf(output_path,"%.255ssmall_turn_left_bank_to_gentle_up%s",output_dir,suffix);
		write_track_section(context,&(track_list.small_turn_left_bank_to_gentle_up),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
		sprintf(output_path,"%.255ssmall_turn_right_bank_to_gentle_up%s",output_dir,suffix);
		write_track_section(context,&(track_list.small_turn_right_bank_to_gentle_up),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	}

	if(groups&TRACK_GROUP_MEDIUM_HALF_LOOPS)
	{
		sprintf(output_path,"%.255smedium_half_loop_left%s",output_dir,suffix);
		write_track_section(context,&(track_list.medium_half_loop_left),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
		sprintf(output_path,"%.255smedium_half_loop_right%s",output_dir,suffix);
		write_track_section(context,&(track_list.medium_half_loop_right),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	}
	if(groups&TRACK_GROUP_LARGE_HALF_LOOPS)
	{
		sprintf(output_path,"%.255slarge_half_loop_left%s",output_dir,suffix);
		write_track_section(context,&(track_list.large_half_loop_left),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
		sprintf(output_path,"%.255slarge_half_loop_right%s",output_dir,suffix);
		write_track_section(context,&(track_list.large_half_loop_right),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	}
	if(groups&TRACK_GROUP_ZERO_G_ROLLS)
	{
		sprintf(output_path,"%.255szero_g_roll_left%s",output_dir,suffix);
		write_track_section(context,&(track_list.zero_g_roll_left),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
		sprintf(output_path,"%.255szero_g_roll_right%s",output_dir,suffix);
		write_track_section(context,&(track_list.zero_g_roll_right),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
		sprintf(output_path,"%.255slarge_zero_g_roll_left%s",output_dir,suffix);
		write_track_section(context,&(track_list.large_zero_g_roll_left),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
		sprintf(output_path,"%.255slarge_zero_g_roll_right%s",output_dir,suffix);
		write_track_section(context,&(track_list.large_zero_g_roll_right),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	}

	if(groups&TRACK_GROUP_SMALL_SLOPE_TRANSITIONS)
	{
		sprintf(output_path,"%.255ssmall_flat_to_steep_up%s",output_dir,suffix);
		write_track_section(context,&(track_list.small_flat_to_steep_up),track_type,base_dir,mask_dir,output_path,sprites,subtype,subtype ==TRACK_SUBTYPE_LIFT ? small_flat_to_steep_up_chain : NULL);
		sprintf(output_path,"%.255ssmall_steep_to_flat_up%s",output_dir,suffix);
		write_track_section(context,&(track_list.small_steep_to_flat_up),track_type,base_dir,mask_dir,output_path,sprites,subtype,subtype ==TRACK_SUBTYPE_LIFT ? small_steep_to_flat_up_chain : NULL);
		sprintf(output_path,"%.255ssmall_flat_to_steep_up_diag%s",output_dir,suffix);
		write_track_section(context,&(track_list.small_flat_to_steep_up_diag),track_type,base_dir,mask_dir,output_path,sprites,subtype,subtype ==TRACK_SUBTYPE_LIFT ? small_flat_to_steep_up_diag_chain : NULL);
		sprintf(output_path,"%.255ssmall_steep_to_flat_up_diag%s",output_dir,suffix);
		write_track_section(context,&(track_list.small_steep_to_flat_up_diag),track_type,base_dir,mask_dir,output_path,sprites,subtype,subtype ==TRACK_SUBTYPE_LIFT ? small_steep_to_flat_up_diag_chain : NULL);
	}

	if(groups&TRACK_GROUP_LARGE_SLOPED_TURNS)
	{
		sprintf(output_path,"%.255slarge_turn_left_to_diag_gentle_up%s",output_dir,suffix);
		write_track_section(context,&(track_list.large_turn_left_to_diag_gentle_up),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
		sprintf(output_path,"%.255slarge_turn_right_to_diag_gentle_up%s",output_dir,suffix);
		write_track_section(context,&(track_list.large_turn_right_to_diag_gentle_up),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
		sprintf(output_path,"%.255slarge_turn_left_to_orthogonal_gentle_up%s",output_dir,suffix);
		write_track_section(context,&(track_list.large_turn_left_to_orthogonal_gentle_up),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
		sprintf(output_path,"%.255slarge_turn_right_to_orthogonal_gentle_up%s",output_dir,suffix);
		write_track_section(context,&(track_list.large_turn_right_to_orthogonal_gentle_up),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	}

	if(groups&TRACK_GROUP_LARGE_BANKED_SLOPED_TURNS)
	{
		sprintf(output_path,"%.255sgentle_up_to_gentle_up_left_bank_diag%s",output_dir,suffix);
		write_track_section(context,&(track_list.gentle_up_to_gentle_up_left_bank_diag),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
		sprintf(output_path,"%.255sgentle_up_to_gentle_up_right_bank_diag%s",output_dir,suffix);
		write_track_section(context,&(track_list.gentle_up_to_gentle_up_right_bank_diag),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
		sprintf(output_path,"%.255sgentle_up_left_bank_to_gentle_up_diag%s",output_dir,suffix);
		write_track_section(context,&(track_list.gentle_up_left_bank_to_gentle_up_diag),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
		sprintf(output_path,"%.255sgentle_up_right_bank_to_gentle_up_diag%s",output_dir,suffix);
		write_track_section(context,&(track_list.gentle_up_right_bank_to_gentle_up_diag),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
		sprintf(output_path,"%.255sleft_bank_to_gentle_up_left_bank_diag%s",output_dir,suffix);
		write_track_section(context,&(track_list.left_bank_to_gentle_up_left_bank_diag),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
		sprintf(output_path,"%.255sright_bank_to_gentle_up_right_bank_diag%s",output_dir,suffix);
		write_track_section(context,&(track_list.right_bank_to_gentle_up_right_bank_diag),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
		sprintf(output_path,"%.255sgentle_up_left_bank_to_left_bank_diag%s",output_dir,suffix);
		write_track_section(context,&(track_list.gentle_up_left_bank_to_left_bank_diag),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
		sprintf(output_path,"%.255sgentle_up_right_bank_to_right_bank_diag%s",output_dir,suffix);
		write_track_section(context,&(track_list.gentle_up_right_bank_to_right_bank_diag),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
		sprintf(output_path,"%.255sgentle_up_left_bank_diag%s",output_dir,suffix);
		write_track_section(context,&(track_list.gentle_up_left_bank_diag),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
		sprintf(output_path,"%.255sgentle_up_right_bank_diag%s",output_dir,suffix);
		write_track_section(context,&(track_list.gentle_up_right_bank_diag),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
		sprintf(output_path,"%.255sflat_to_gentle_up_left_bank_diag%s",output_dir,suffix);
		write_track_section(context,&(track_list.flat_to_gentle_up_left_bank_diag),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
		sprintf(output_path,"%.255sflat_to_gentle_up_right_bank_diag%s",output_dir,suffix);
		write_track_section(context,&(track_list.flat_to_gentle_up_right_bank_diag),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
		sprintf(output_path,"%.255sgentle_up_left_bank_to_flat_diag%s",output_dir,suffix);
		write_track_section(context,&(track_list.gentle_up_left_bank_to_flat_diag),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
		sprintf(output_path,"%.255sgentle_up_right_bank_to_flat_diag%s",output_dir,suffix);
		write_track_section(context,&(track_list.gentle_up_right_bank_to_flat_diag),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
		sprintf(output_path,"%.255slarge_turn_left_bank_to_diag_gentle_up%s",output_dir,suffix);
		write_track_section(context,&(track_list.large_turn_left_bank_to_diag_gentle_up),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
		sprintf(output_path,"%.255slarge_turn_right_bank_to_diag_gentle_up%s",output_dir,suffix);
		write_track_section(context,&(track_list.large_turn_right_bank_to_diag_gentle_up),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
		sprintf(output_path,"%.255slarge_turn_left_bank_to_orthogonal_gentle_up%s",output_dir,suffix);
		write_track_section(context,&(track_list.large_turn_left_bank_to_orthogonal_gentle_up),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
		sprintf(output_path,"%.255slarge_turn_right_bank_to_orthogonal_gentle_up%s",output_dir,suffix);
		write_track_section(context,&(track_list.large_turn_right_bank_to_orthogonal_gentle_up),track_type,base_dir,mask_dir,output_path,sprites,subtype,NULL);
	}
	return 0;
}

int write_track_type(context_t* context,track_type_t* track_type,json_t* sprites,const char* base_dir, const char* mask_dir,const char* output_dir)
{
	track_list_t track_list=track_list_default;

	std::string mask_type_dir = std::format("{}{}/", mask_dir, track_type->masks_name);

	write_track_subtype(context,track_type,track_list,sprites,base_dir,mask_type_dir.c_str(),output_dir,TRACK_SUBTYPE_DEFAULT);
	if(track_type->flags&TRACK_HAS_LIFT)
	{
		write_track_subtype(context,track_type,track_list,sprites,base_dir,mask_type_dir.c_str(),output_dir,TRACK_SUBTYPE_LIFT);
	}
	return 0;
}
